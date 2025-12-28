function [k,v_] = solveBlochFluid(w,k_target,Re,YPSub,a,b,shape,Nx,Ny,Nm,mesh_type)
%SOLVEBLOCHFLUID Solves spatial Orr-Sommerfeld problem for flow over a periodic surface.
%
%   Computes the Bloch wavenumbers and mode shapes for a fluid flowing over
%   a periodic surface (defined by admittance YPSub and unit cell shape).
%   The method uses a hybrid Fourier (streamwise/spanwise) + FEM (wall-normal) 
%   discretization.
%
%   Inputs:
%       w         : Non-dimensional frequency (omega)
%       k_target  : Target non-dimensional wavenumber (guess for eigenvalue solver)
%       Re        : Reynolds number
%       YPSub     : Non-dimensional admittance of the bottom wall
%       a         : Unit cell dimensions (periodicity size)
%       b         : Square PSub side length 
%       shape     : 'none', 'fullspan', 'square', or 'hexagon'
%       Nx        : Number of Fourier modes (streamwise/spanwise resolution)
%       Ny        : Number of Finite Elements (wall-normal resolution)
%       Nm        : Number of Bloch modes to compute
%       mesh_type : (Optional) Grid stretching: 'uniform', 'cosine' (default)
%
%   Outputs:
%       k         : Vector of computed Bloch wavenumbers (eigenvalues)
%       v_        : Computed Bloch mode shapes (eigenvectors)

%% 1. Validate inputs and set defaults
if ~exist('mesh_type','var') || isempty(mesh_type)
    mesh_type = 'cosine'; 
end
valid_shape = {'none','fullspan','square','hexagon'};
valid_mesh_type = {'uniform','cosine'};
shape = validatestring(shape, valid_shape);
mesh_type = validatestring(mesh_type, valid_mesh_type);

%% 2. Discretization setup

% Compute FEM matrices for wall-normal direction (normalized between y=0 and y=2)
[y,Nf,Bf,Hf,Cf,~,C2f,Df,~] = computeFEMatrices(Ny,mesh_type);
n_dof_fem = size(Nf,1);

% Get Fourier coefficients for the unit cell step function
% H_ : Step function coefficients, a_ : Fourier wavenumbers
[H_,a_] = getStepFourier(a,b,shape,Nx,true);
n_fourier = size(a_,1);

%% 3. Matrices coefficients calculation

% Get fluid matrices coefficients
[C_] = getFluidMatricesCoefficients(a_,Re);

% Get boundary conditions coefficients
[Bw_,Bk_,BY_] = getBoundaryCoefficients(H_,a_,y,Re);

%% 4. Matrices assembly

% Generate connectivities matrices for the global system
[Tn_,P] = getFluidConnectivities(n_dof_fem,n_fourier);

% Compute fluid matrices
fem_matrices = {Nf,Cf,C2f,Bf,Df,Hf};
A_ = computeFluidMatrices(Tn_,fem_matrices,C_,Bw_,Bk_,BY_,P,w,'spatial',YPSub);

%% 5. Eigenvalues problem resolution

% Solve for eigenvalues close to k_target
[V,alpha] = neigs(A_,Nm,k_target);

% Assign outputs
k = alpha;
v_ = reshape(P*V,size(Tn_,1),size(Tn_,2),size(V,2));

end