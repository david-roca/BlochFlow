function [dof_conn,con_mat] = getFluidConnectivities(dof_mode,N_modes)
%GETFLUIDCONNECTIVITIES Generates DOF maps and applies boundary conditions.
%
%   Constructs the global connectivity map and the constraint projection matrix.
%   Assumes the top boundary (last node) is a rigid wall (v=0, dv/dy=0).
%   The bottom boundary (node 1) is left free for the FSI interface.
%
%   Inputs:
%       dof_mode : Number of vertical DOFs (nodes*2).
%       N_modes  : Number of Fourier modes being simulated.
%
%   Outputs:
%       dof_conn : [dof_mode x N_modes] Global index map.
%       con_mat  : [N_total x N_reduced] Projection matrix used to 
%                  reduce the system matrices (removes constrained DOFs).
%                  x_full = P * x_reduced

%% 1. Global DOF map

% Total system size
N_dofs = dof_mode * N_modes;

% Reshape linear index 1:N into a matrix [wall-normal dof, mode index]
dof_conn = reshape(1:N_dofs, dof_mode, N_modes);

%% 2. Identify constrained DOFs (top wall)

% Last two DOFs (value and slope at the top wall) constrained for every mode 
ind_con = dof_conn(end-1:end, :);

%% 3. Build projection matrix imposing constraints

% Find the indices that remain free (set difference)
ind_free = setdiff(1:N_dofs, ind_con(:));

% Number of reduced DOFs
n_red = length(ind_free);

% Construct sparse matrix directly 
con_mat = sparse(ind_free, 1:n_red, 1, N_dofs, n_red);
    
end