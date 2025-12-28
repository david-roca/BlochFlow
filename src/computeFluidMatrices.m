function [A_] = computeFluidMatrices(dof_conn,FE_mat,C_,Bw_,Bk_,BY_,P_con,param,method,Y)
%COMPUTEFLUIDMATRICES Assembles the global stability system matrices.
%
%   Constructs the matrices A_i for the eigenvalue problem:
%   (A_1 + lambda*A_2 + lambda^2*A_3 + ...) * x = 0
%
%   Inputs:
%       dof_map : [N_dof x N_modes] Connectivity map.
%       FE_mats : Cell array of FE integrals (Mass, Stiffness, etc.).
%                 {1}:N, {2}:C, {3}:C2, {4}:B, {5}:D, {6}:H
%       C_      : Physical coefficients from getFluidMatricesCoefficients.
%       Bw_     : Boundary condition matrices (frequency terms).
%       Bk_     : Boundary condition matrices (wavenumber terms).
%       BY_     : Boundary condition matrices (admittance terms).
%       P_con   : Projection matrix for top-wall constraints.
%       param   : The fixed parameter:
%                 - Wavenumber 'k' if method is 'temporal'.
%                 - Frequency 'w' if method is 'spatial'.
%       method  : 'temporal' or 'spatial'.
%       Y       : Admittance value (scalar) for the compliant wall.
%
%   Outputs:
%       A_      : Cell array of global sparse matrices.

%% 1. Dimensions

% Get dimensions
nc = length(C_);
nm = length(FE_mat);
nf = size(dof_conn,2);
ni = size(dof_conn,1);
Ndof = nf*ni;

%% 2. Compute scalar weights

switch method
    case 'temporal' % Solve for frequency given wavenumber A0 + w*A1 = 0
        nd = 2; % Polynomial degree + 1 (number of matrices)
        phi = zeros(nf,nm,nd);
        for i = 1:nc
            phi(:,i,1) = sum(C_{i}(:,[1,2,4,6,8]).*[1,param,param^2,param^3,param^4],2);
            phi(:,i,2) = sum(C_{i}(:,[3,5,7]).*[1,param,param^2],2);
        end
        phi(:,nm,1) = 1;
    case 'spatial' % Solve for wavenumber given frequency A0 + ... + k^4*A4 = 0
        nd = 5; % Polynomial degree + 1 (number of matrices)
        phi = zeros(nf,nm,nd);
        for i = 1:length(C_)
            phi(:,i,1) = sum(C_{i}(:,[1,3]).*[1,param],2);
            phi(:,i,2) = sum(C_{i}(:,[2,5]).*[1,param],2);
            phi(:,i,3) = sum(C_{i}(:,[4,7]).*[1,param],2);
            phi(:,i,4) = C_{i}(:,6);
            phi(:,i,5) = C_{i}(:,8);
        end
        phi(:,nm,1) = 1;
    otherwise
        error('Method must be "temporal" or "spatial"');
end

%% 3. Bulk matrix assembly

A_ = cell(1,nd);        
for i = 1:nd
    % Initialize sparse matrix
    A_{i} = sparse(Ndof,Ndof);
    for n = 1:nf
        % Sum weighted FE matrices
        idx = dof_conn(:,n);
        A_loc = sparse(ni,ni);
        for j = 1:nm
            A_loc = A_loc + phi(n,j,i)*FE_mat{j};
        end
        % Add to global system
        A_{i}(idx,idx) = A_{i}(idx,idx) + A_loc;
    end
end

%% 4. Apply boundary conditions matrices

% a) Clear first two rows for all modes (where BCC equations will be placed) 
for i = 1:length(A_)
    A_{i}(dof_conn(1:2,:),:) = 0;
end

% b) Extend A_ if BCs require higher polynomial order
for i = (length(A_)+1):length(Bk_)
    A_{i} = sparse(Ndof,Ndof);
end

% c) Assemble BC matrices
idof = repmat(repmat(dof_conn(1:2,:),4,1),1,nf);
jdof = repelem(repelem(dof_conn(1:4,:),2,1),1,nf);
switch method
    case 'temporal'
        % w is the polynomial variable (frequency)
        for i = 1:length(Bw_)
            A_{i} = A_{i} + sparse(idof(:),jdof(:),Bw_{i}(:),Ndof,Ndof);
        end
        % k = param (wavenumber) - that is why it only adds to matrix A_1        
        for i = 1:length(Bk_) 
            A_{1} = A_{1} + sparse(idof(:),jdof(:),Bk_{i}(:)*param^(i-1),Ndof,Ndof);
        end
        for i = 1:length(BY_)
            A_{1} = A_{1} + sparse(idof(:),jdof(:),BY_{i}(:)*(param^(i-1))*Y,Ndof,Ndof);
        end
    case 'spatial'
        % w = param (frequency) - that is why it only adds to matrix A_1 
        for i = 1:length(Bw_)
            A_{1} = A_{1} + sparse(idof(:),jdof(:),Bw_{i}(:)*param^(i-1),Ndof,Ndof);
        end
        % k is the polynomial variable (wavenumber)
        for i = 1:length(Bk_)
            A_{i} = A_{i} + sparse(idof(:),jdof(:),Bk_{i}(:),Ndof,Ndof);
        end
        for i = 1:length(BY_)
            A_{i} = A_{i} + sparse(idof(:),jdof(:),BY_{i}(:)*Y,Ndof,Ndof);
        end
end

%% 5. Apply prescribed DOFs

% Remove constrained DOFs rows and columns
for i = 1:length(A_)
    A_{i} = P_con'*A_{i}*P_con;
end

end