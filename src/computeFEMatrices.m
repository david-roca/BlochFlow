function [y,Nf,Bf,Hf,Cf,C1f,C2f,Df,Ef] = computeFEMatrices(N_elem,mesh_type)
%COMPUTEFEMATRICES Assembles 1D FEM matrices for the Orr-Sommerfeld equation.
%
%   Computes the system matrices using Hermite Cubic elements (C1 continuity).
%   The domain is normalized to y in [0, 2].
%
%   Inputs:
%       N_elem    : Number of elements.
%       mesh_type : 'uniform' (linear) or 'cosine' (Chebyshev-like clustering).
%
%   Outputs:
%       y         : [N_nodes x 1] Nodal coordinates.
%       Nf        : "Mass" matrix (Integral of N*N).
%       Bf        : "Stiffness" matrix (Integral of N'*N').
%       Hf        : "Biharmonic" matrix (Integral of N''*N'').
%       Cf        : Convection matrix 0 (Base flow U * Mass).
%       C1f       : Convection matrix 1 (Base flow shear U' * Mass).
%       C2f       : Convection matrix 2 (Base flow curvature U'' * Mass).
%       Df        : Mixed matrix (Integral of U * N * N'').
%       Ef        : Mixed matrix (Integral of N * N').
%
%   Dependencies:
%       - getShapeFun(xi) : Returns Hermite shape functions and derivatives.
%       - getBaseFlow(y)  : Returns base flow profile (U, U', U'').

%% 1. Mesh generation

% Check inputs
if ~exist('mesh_type','var') || isempty(mesh_type)
    mesh_type = 'cosine'; 
end

% Discretization
switch lower(mesh_type)
    case 'uniform'
        y = linspace(0,2,N_elem+1);
    case 'cosine'
        y = 1-cos(pi*linspace(0,1,N_elem+1));
    otherwise
        error('Unknown mesh type. Use "uniform" or "cosine".');
end
nod_conn = [1:N_elem;2:N_elem+1]; % Nodal connectivities
nod_elem = size(nod_conn,1);      % Nodes per element
dof_node = 2;                     % DOF 1: stream function, DOF 2: derivative

% Gauss points data
ngauss = 5;
[xg,wg] = getGaussPoints(ngauss);

% Build element matrices
dims = [nod_elem*dof_node,nod_elem*dof_node,N_elem];
Ne = zeros(dims);
Be = zeros(dims);
He = zeros(dims);
Ce = zeros(dims);
C1e = zeros(dims);
C2e = zeros(dims);
De = zeros(dims);
Ee = zeros(dims);

%% Element assembly loop
for i = 1:ngauss
    % Get Hermite shape functions at Gauss point i
    [N_val,N_der,N_dder] = getShapeFun(xg(i));
    
    % Compute Jacobian
    y_elem = y(nod_conn);
    J = 1./(1-sum(N_der([2,4]))).*N_der([1,3])'*y_elem;
    detJ = J;
    invJ = 1./J;

    % Scaling operator from physical slopes to element coordinates slopes
    Nop = [ones(1,N_elem);J;ones(1,N_elem);J];
    
    % Operators to transform derivatives
    Bop = invJ.*Nop;
    Hop = invJ.*Bop;

    % Transform shape functions to physical space
    N_phy = Nop.*N_val;
    B_phy = Bop.*N_der;
    H_phy = Hop.*N_dder;

    % Get base flow properties at physical Gauss points
    yg = sum([y_elem(1,:);ones(1,N_elem);y_elem(2,:);ones(1,N_elem)].*N_phy,1); % Retrieve Gauss points coordinates
    [U_,Up_,Upp_] = getBaseFlow(yg);

    % Permuted terms
    dy = wg(i)*permute(detJ,[1,3,2]);
    N_T = permute(N_phy, [1,3,2]);
    N_  = permute(N_phy, [3,1,2]);
    B_T = permute(B_phy, [1,3,2]);
    B_  = permute(B_phy, [3,1,2]);
    H_T = permute(H_phy, [1,3,2]);
    H_  = permute(H_phy, [3,1,2]);
    
    % Calculate element matrices
    Ne  = Ne  + dy.*N_T.*N_;
    Be  = Be  + dy.*B_T.*B_;
    He  = He  + dy.*H_T.*H_;
    Ce  = Ce  + dy.*permute(U_,  [1,3,2]).*N_T.*N_;
    C1e = C1e + dy.*permute(Up_, [1,3,2]).*N_T.*N_;
    C2e = C2e + dy.*permute(Upp_,[1,3,2]).*N_T.*N_;
    De  = De  + dy.*permute(U_,  [1,3,2]).*N_T.*H_;
    Ee  = Ee  + dy.*N_T.*B_;
end

%% 4. Global assembly

% Total number of DOFs
ndof = dof_node*(N_elem+1);

% Generate DOF indices for sparse matrices assembly
dof_conn = reshape(dof_node*(permute(nod_conn,[3,1,2])-1)+(1:dof_node)', ...
                   dof_node*nod_elem,N_elem);
idof = repmat(dof_conn,dof_node*nod_elem,1); 
jdof = repelem(dof_conn,dof_node*nod_elem,1);

% Sparse matrices
Nf  = sparse(idof(:),jdof(:),Ne(:) ,ndof,ndof);  
Bf  = sparse(idof(:),jdof(:),Be(:) ,ndof,ndof);  
Hf  = sparse(idof(:),jdof(:),He(:) ,ndof,ndof);  
Cf  = sparse(idof(:),jdof(:),Ce(:) ,ndof,ndof);  
C1f = sparse(idof(:),jdof(:),C1e(:),ndof,ndof);  
C2f = sparse(idof(:),jdof(:),C2e(:),ndof,ndof);  
Df  = sparse(idof(:),jdof(:),De(:) ,ndof,ndof);  
Ef  = sparse(idof(:),jdof(:),Ee(:) ,ndof,ndof); 

end

%--------------------------------------------------------------------------
% LOCAL HELPER FUNCTIONS
%--------------------------------------------------------------------------
function [xg,wg] = getGaussPoints(ngauss)
%GETGAUSSPOINTS Returns coordinates and weights for Gauss-Legendre quadrature.
    xg = zeros(1,ngauss);
    wg = zeros(1,ngauss);
    switch ngauss
        case 1
            xg(1) = 0;
            wg(1) = 2;
        case 2
            xg([1,2]) = [-1,1]/sqrt(3);
            wg([1,2]) = [1,1];
        case 3
            xg([1,3]) = [-1,1]*sqrt(3/5);
            wg([1,3]) = [1,1]*5/9;
            wg(2) = 8/9;
        case 4
            xg([1,4]) = [-1,1]*sqrt(3/7+2/7*sqrt(6/5));
            xg([2,3]) = [-1,1]*sqrt(3/7-2/7*sqrt(6/5));
            wg([1,4]) = [1,1]*(18-sqrt(30))/36;
            wg([2,3]) = [1,1]*(18+sqrt(30))/36;
        case 5
            xg([1,5]) = [-1,1]*sqrt(5+2*sqrt(10/7))/3;
            xg([2,4]) = [-1,1]*sqrt(5-2*sqrt(10/7))/3;
            wg([1,5]) = [1,1]*(322-13*sqrt(70))/900;
            wg([2,4]) = [1,1]*(322+13*sqrt(70))/900;
            wg(3) = 128/225;
        case 6
            xg([1,6]) = [-1,1]*0.9324695142031521;
            xg([2,5]) = [-1,1]*0.6612093864662645;
            xg([3,4]) = [-1,1]*0.2386191860831969;
            wg([1,6]) = [1,1]*0.1713244923791704;
            wg([2,5]) = [1,1]*0.3607615730481386;
            wg([3,4]) = [1,1]*0.4679139345726910;
        case 7
            xg([1,7]) = [-1,1]*0.9491079123427585;
            xg([2,6]) = [-1,1]*0.7415311855993945;
            xg([3,5]) = [-1,1]*0.4058451513773972;
            wg([1,7]) = [1,1]*0.1294849661688697;
            wg([2,6]) = [1,1]*0.2797053914892766;
            wg([3,5]) = [1,1]*0.3818300505051189;
            wg(4) = 0.4179591836734694;
        case 8
            xg([1,8]) = [-1,1]*0.9602898564975363;
            xg([2,7]) = [-1,1]*0.7966664774136267;
            xg([3,6]) = [-1,1]*0.5255324099163290;
            xg([4,5]) = [-1,1]*0.1834346424956498;
            wg([1,8]) = [1,1]*0.1012285362903763;
            wg([2,7]) = [1,1]*0.2223810344533745;
            wg([3,6]) = [1,1]*0.3137066458778873;
            wg([4,5]) = [1,1]*0.3626837833783620;
        otherwise
            error('Gauss points > 5 not implemented in this cleaned version (but can be added).');
    end
end

