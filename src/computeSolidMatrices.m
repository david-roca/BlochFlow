function [K,M,C,K_uc,M_uc,C_uc] = computeSolidMatrices(l_cell,A,rho,E,q_damp,m_res,k_res,N_cell,Ny)
%COMPUTESOLIDMATRICES Assembles FE matrices for a 1D rod with optional resonators.
%
%   Constructs the Stiffness (K), Mass (M), and Damping (C) matrices for a
%   periodic rod structure. The structure can be a simple Phononic Crystal (PnC)
%   or a Metamaterial (MM) with internal resonators.
%
%   Inputs:
%       l_cell  : Length of the unit cell (scalar or [l1, l2])
%       A       : Cross-section area (scalar or [A1, A2])
%       rho     : Material density (scalar or [rho1, rho2])
%       E       : Young's modulus (scalar or [E1, E2])
%       q_damp  : Rayleigh damping coefficients [alpha, beta] (C = alpha*M + beta*K)
%       m_res   : Resonator mass (scalar or array of length N_cells). Empty [] for PnC.
%       k_res   : Resonator stiffness (scalar or array of length N_cells).
%       N_cell  : Number of unit cells.
%       Ny      : (Optional) Number of Finite Elements per segment. Default: 1.
%
%   Outputs:
%       K, M, C       : Global sparse system matrices (with BC applied)
%       K_uc, M_uc... : Matrices for a single unit cell (useful for Bloch analysis)

%% 1. Validate inputs and set defaults

if ~exist('Ny','var') || isempty(Ny), Ny = 1; end
if isscalar(l_cell), l_cell = [l_cell/2,l_cell/2]; end
if isscalar(A), A = [A,A]; end
if isscalar(rho), rho = [rho,rho]; end
if isscalar(E), E = [E,E]; end

%% 2. Define connectivities

% Determine whether we have resonators or not
resonators = ~isempty(m_res);

if ~resonators
    % --- PnC (no internal resonators) ---
    % Nodes for segments 1 and 2
    ind1 = (1:(Ny+1))'        + (0:((2*Ny)):((2*Ny)*(N_cell-1)));
    ind2 = ((Ny+1):(2*Ny+1))' + (0:((2*Ny)):((2*Ny)*(N_cell-1)));
    indr = zeros(0,N_cell); % No resonator nodes
else 
    % --- Metamaterial (with resonators) ---
    % Expand resonator properties if scalar
    if isscalar(m_res)
        if ~exist('N_cell','var') || isempty(N_cell), N_cell = 1; end
        m_res = m_res*ones(1,N_cell);
        k_res = k_res*ones(1,N_cell);
    else
        if ~exist('N_cell','var') || isempty(N_cell)
            N_cell = length(m_res);
        elseif N_cell ~= length(m_res)
            error('Number of unit cells provided does not match the dimensions of the internal resonators'' mass array');
        end
    end
    % Nodes for segments 1 and 2 and resonators
    ind1 = (1:(Ny+1))'             + (0:(2*Ny+1):((2*Ny+1)*(N_cell-1)));
    ind2 = [Ny+1,(Ny+3):(2*Ny+2)]' + (0:(2*Ny+1):((2*Ny+1)*(N_cell-1)));
    indr = (Ny+2)                  + (0:(2*Ny+1):((2*Ny+1)*(N_cell-1)));
end

%% 3. Element matrices computation

% Calculate number of degrees of freedom (DOFs) per unit cell block
dof_unit_cell = size(ind1,1)+size(ind2,1)+size(indr,1)-1;

% Preallocate element matrices in 3D array
K_ = zeros(dof_unit_cell,dof_unit_cell,N_cell);
M_ = zeros(dof_unit_cell,dof_unit_cell,N_cell);

% Standard 1D rod element matrices
K_rod = @(E_val, A_val, L_val) (E_val*A_val/L_val)*[1, -1; -1, 1];
M_rod = @(rho_val, A_val, L_val) (rho_val*A_val*L_val/6)*[2, 1; 1, 2];

% Loop over unit cells 
for i = 1:N_cell
    % Segment 1 assembly
    for j = 1:size(ind1,1)-1
        idx = ind1([j,j+1],1); % Local indices
        l_ele = l_cell(1)/Ny;
        K_(idx,idx,i) = K_(idx,idx,i) + K_rod(E(1), A(1), l_ele);
        M_(idx,idx,i) = M_(idx,idx,i) + M_rod(rho(1), A(1), l_ele);
    end
    % Segment 2 assembly
    for j = 1:size(ind2,1)-1
        idx = ind2([j,j+1],1);
        l_ele = l_cell(2)/Ny;
        K_(idx,idx,i) = K_(idx,idx,i) + K_rod(E(2), A(2), l_ele);
        M_(idx,idx,i) = M_(idx,idx,i) + M_rod(rho(2), A(2), l_ele);
    end
    % Resonator assembly
    if resonators
        idx = [ind1(end,1),indr(1)]; % Resonator attached to mid-node
        K_(idx,idx,i) = K_(idx,idx,i) + k_res(i)*[1, -1; -1, 1];
        M_(idx,idx,i) = M_(idx,idx,i) + m_res(i)*[0, 0; 0, 1];
    end
end

% Rayleigh damping matrix
C_ = q_damp(1)*M_+q_damp(2)*K_;

%% 4. Unit cell outputs

% Retrieve element matrices
if nargout > 3
    K_uc = K_(:,:,1);
    M_uc = M_(:,:,1);
    C_uc = C_(:,:,1);
end

%% 5. Global assembly

% Total number of DOFs
ndof = (dof_unit_cell-1)*N_cell+1;

% Indices mapping vectors
ind_vec = reshape([ind1;indr;ind2(2:end,:)],[],1);
Idof = repmat(ind_vec,dof_unit_cell,1);
Jdof = repelem(ind_vec,dof_unit_cell,1);

% Sparse matrices assembly
K = sparse(Idof(:),Jdof(:),K_(:),ndof,ndof);
M = sparse(Idof(:),Jdof(:),M_(:),ndof,ndof);
C = sparse(Idof(:),Jdof(:),C_(:),ndof,ndof);

%% 6. Apply boundary conditions

% Remove last row/col (fixed end)
K = K(1:end-1,1:end-1);
M = M(1:end-1,1:end-1);
C = C(1:end-1,1:end-1);

end