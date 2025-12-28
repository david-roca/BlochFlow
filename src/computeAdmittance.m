function Y = computeAdmittance(w,K,M,C,A)
%COMPUTEADMITTANCE Computes the wall admittance (Y = v/p) for the solid domain.
%
%   Y = COMPUTEADMITTANCE(w, K, M, C, A) calculates the acoustic admittance
%   of the solid surface by reducing the full system matrices to the
%   fluid-structure interface node via static condensation (Schur complement).
%
%   Inputs:
%       w        : [N x 1] Array of angular frequencies (rad/s).
%       K, M, C  : Global Stiffness, Mass, and Damping matrices.
%       Area     : Cross-sectional area of the fluid-solid interface.
%
%   Outputs:
%       Y        : [N x 1] Complex admittance values (velocity / pressure).

% Get number of DOFs in solid domain

%% 1. Initialization
N = size(K,1);
ind1 = 1;   % First node is the interface with the flow
indi = 2:N; % Internal nodes
Y = zeros(size(w));

%% 2. Loop through frequencies
for i = 1:length(w)
    
    % Compute dynamic stiffness matrix
    D = K-1i*w(i)*C-w(i)^2*M;

    % Schur complement condensation
    D_eff = D(ind1,ind1)-D(ind1,indi)*(D(indi,indi)\D(indi,ind1));

    % Compute admittance Y = v/p
    % 1) force is F = -p*A (negative because p pushes IN and u is positive OUT)
    % 2) displacement is u = F/D_eff
    % 3) velocity is v = -i*w*u [time harmonic exp(-i*w*t)]
    % 4) admittance is Y = v/p = (i*w*A)/D_eff
    Y(i) = (1i*w(i)*A)./D_eff;

end
end