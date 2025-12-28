function [Bw_,Bk_,BY_] = getBoundaryCoefficients(H_,k_,y,Re)
%GETBOUNDARYCOEFFICIENTS Assembles FSI boundary condition matrices.
%
%   Constructs the matrices representing the linearized boundary conditions
%   at the fluid-solid interface (y=0).
%   
%   The function computes the coefficients for the Polynomial Eigenvalue 
%   Problem (PEP) associated with the pressure and velocity coupling.
%
%   Inputs:
%       H_   : [N x 1] Fourier coefficients of the wall step function.
%       k_   : [N x 2] Fourier wavenumbers [alpha, beta] for each mode.
%       y    : Nodal coordinates (used for Jacobian calculation).
%       Re   : Reynolds number.
%
%   Outputs:
%       Bw_  : Cell array of matrices associated with frequency terms.
%       Bk_  : Cell array of matrices associated with wavenumber terms.
%       BY_  : Cell array of matrices associated with admittance terms.
%
%   The output matrices are 4D arrays: [BC_row, elem_DOF, row_mode, col_mode]
%       - BC_row 1: Dynamic condition (pressure)
%       - BC_row 2: Kinematic condition (velocity)

%% 1. Pre-computations

% Generate convolution coefficients representing the operator expansion
[F,F_] = getConvolutionCoefficients(k_);

% Jacobian at the bottom wall
J = (y(2)-y(1))/2;

% Shape function third derivative at wall (needed for pressure gradient)
[~,~,~,N_ddder] = getShapeFun(-1);

% Transform derivative to physical space: d^3/dy^3 = (1/J)^3*d^3/dx^3
N3_phy = (1/J^3)*N_ddder.*[1;J;1;J];

% Coupling interaction matrix
c_ = 2*H_*(H_'/H_(1));

%% 2. Initialization of cell arrays

% Initialize matrices
N = size(k_,1);
N_pol = length(F);
Bw_ = cell(2,1);
Bk_ = cell(N_pol,1);
BY_ = cell(N_pol,1);
dim = [2,4,N,N];
for i = 1:length(Bw_), Bw_{i} = zeros(dim); end
for i = 1:length(Bk_), Bk_{i} = zeros(dim); end
for i = 1:length(BY_), BY_{i} = zeros(dim); end

% Get base flow shear at wall (Up_wall = dU/dy)
[~,Up_wall] = getBaseFlow(0); % Base flow velocity derivative at bottom wall

%% 3. Populate matrices

% Populate matrices
for n = 1:N
    for m = 1:N
        if n==m
            % Kinematic condition (v_fluid = v_solid)
            Bw_{2}(:,:,n,m) = [
                0, 0, 0, 0;
                0, 1, 0, 0;
            ];
            Bk_{2}(:,:,n,m) = Up_wall*[
                0, 0, 0, 0;
                1, 0, 0, 0;
            ];
            Bk_{1}(:,:,n,m) = (k_(n,1)*Up_wall)*[
                0, 0, 0, 0;
                1, 0, 0, 0;
            ];
        end
        % Dynamic condition (pressure balance)
        for i = 1:N_pol
            if n==m
                Bk_{i}(:,:,n,m) = Bk_{i}(:,:,n,m) + F(i)*[
                    1, 0, 0, 0;
                    0, 0, 0, 0;
                ];
            end
            BY_{i}(:,:,n,m) = BY_{i}(:,:,n,m) + (F(i)*c_(n,m)/Re)*[
                0, 1, 0, 0;
                0, 0, 0, 0;
            ];
        end
        for i = 1:length(F_{m})
            BY_{i}(:,:,n,m) = BY_{i}(:,:,n,m) - (F_{m}(i)*c_(n,m)/Re)*[
                N3_phy(1), N3_phy(2), N3_phy(3), N3_phy(4);
                        0,         0,         0,         0;
            ];
        end
    end
end

end

% -------------------------------------------------------------------------
% LOCAL HELPER: CONVOLUTION COEFFICIENTS
% -------------------------------------------------------------------------
function [F,F_] = getConvolutionCoefficients(k_)
%GETCONVOLUTIONCOEFFICIENTS Computes characteristic polynomial coefficients.
%
%   Constructs polynomials corresponding to the product of operators:
%   Product_m [ (d/dy)^2 - (alpha_m^2 + beta_m^2) ]
%   This arises from solving the pressure Poisson equation in Fourier space.

% Get the number of Fourier modes
N = size(k_,1);

% Initialize polynomials
F = 1;
F_ = num2cell(ones(N,1));

% Build polynomials through convolution
for m = 1:N
    alpha = k_(m,1);
    beta = k_(m,2);
    % The operator kernel for mode m: (s + alpha)^2 + beta^2
    kernel = [1, 2*alpha, alpha^2 + beta^2];
    % Multiply into global polynomial
    F = conv(F,kernel);
    % Multiply into all other more polynomials
    for n = 1:N
        if n~=m
            F_{n} = conv(F_{n},kernel);
        end
    end
end

% Normalize and flip coefficients order
Fmax = max(abs(F));
F = flip(F)/Fmax;
for i = 1:length(F_)
    F_{i} = flip(F_{i})/Fmax;
end

end