function [C_] = getFluidMatricesCoefficients(k_,Re)
%GETFLUIDMATRICESCOEFFICIENTS Computes scalar coefficients for system assembly.
%
%   Generates the coefficients c_ij that multiply the Finite Element
%   matrices (Mass, Stiffness, Convection, etc.) to form the full 
%   Orr-Sommerfeld/Squire operator.
%
%   Inputs:
%       k_ : [N x 2] Matrix of wavenumbers [alpha, beta].
%       Re : Reynolds number.
%
%   Outputs:
%       C_ : {5 x 1} Cell array of [N x 8] coefficient matrices.
%            Indices correspond to different FE matrices:
%            1: Biharmonic / Viscous terms
%            2: Convection (U) terms
%            3: Shear (U') terms
%            4: Laplacian terms
%            5: Gradient terms

%% 1. Pre-computations

% Retrieve wavenumbers
alpha = k_(:,1);
beta = k_(:,2);

% Squared wavenumber magnitude
k2 = alpha.^2+beta.^2;

% Imaginary Reynolds
iRe = 1i*Re;

% Helper variables
N = length(alpha);
vec1 = ones(N,1);
vec0 = zeros(N,1);

% Compute coefficients
C_{1} = [
    k2.^2,...                   % -
    4*alpha.*k2,...             % k
    -iRe*k2,...                 % w
    6*alpha.^2 + 2*beta.^2,...  % k^2
    -iRe*2*alpha,...            % k*w
    4*alpha,...                 % k^3
    -iRe*vec1,...               % k^2*w
    vec1;                       % k^4
];
C_{2} = [
    iRe*alpha.*k2,...               %
    iRe*(3*alpha.^2 + beta.^2),...  % k
    vec0,...                        % w
    iRe*3*alpha,...                 % k^2
    vec0,...                        % k*w
    iRe*vec1,...                    % k^3
    vec0,...                        % k^2*w
    vec0;                           % k^4
];
C_{3} = [
    iRe*alpha,...  %
    iRe*vec1,...   % k
    vec0,...       % w
    vec0,...       % k^2
    vec0,...       % k*w
    vec0,...       % k^3
    vec0,...       % k^2*w
    vec0;          % k^4
];
C_{4} = [
    2*k2,...       %
    4*alpha,...    % k
    -iRe*vec1,...  % w
    2*vec1,...     % k^2
    vec0,...       % k*w
    vec0,...       % k^3
    vec0,...       % k^2*w
    vec0;          % k^4
];
C_{5} = [
    -iRe*alpha,... %
    -iRe*vec1,...  % k
    vec0,...       % w
    vec0,...       % k^2
    vec0,...       % k*w
    vec0,...       % k^3
    vec0,...       % k^2*w
    vec0;          % k^4
];
end