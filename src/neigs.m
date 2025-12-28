function [V,D] = neigs(A_,k,sigma)
%NEIGS Solves the Polynomial Eigenvalue Problem via linearization.
%
%   Solving: (A_1 + lambda*A_2 + lambda^2*A_3 + ... + lambda^N*A_N) * x = 0
%
%   The problem is converted into a companion matrix form (Generalized 
%   Eigenvalue Problem) of size (N-1)*size(A).
%
%   Inputs:
%       A_     : Cell array {A1, A2, ..., AN} of square matrices.
%       k      : Number of eigenvalues to compute.
%       sigma  : Target shift (e.g., 0, 'sm', or a complex number).
%
%   Outputs:
%       V      : [size(A) x k] Matrix of eigenvectors (physical DOFs only).
%       D      : [k x 1] Vector of sorted eigenvalues.

%% 1. Problem setup

% Degrees of freedom
ndof = size(A_{1},1);

% Polynomial degree
nd = length(A_);

% Size of linearized system
Ndof = (nd-1)*ndof;

%% 2. Build companion matrices

% We solve the system A*Y = lambda*B*Y, with Y = [x; lambda*x; lambda^2*x ...]
A = sparse(Ndof,Ndof);
B = sparse(Ndof,Ndof);

% Row block 1 (physics equation)
id_1 = 1:ndof;
A(id_1,id_1) = A_{1};
for j = 1:nd-1
    id_j = (j-1)*ndof + (1:ndof);
    B(id_1,id_j) = -A_{j+1};
end

% Remaining row blocks (identity shifts) Y_2 = lambda*Y_1, Y_3 = lambda*Y_2, etc.
for j = 2:nd-1
    id_j = (j-1)*ndof + (1:ndof);
    id_j_prev = (j-2)*ndof + (1:ndof);
    A(id_j, id_j) = speye(ndof,ndof);
    B(id_j, id_j_prev) = speye(ndof,ndof);
end

%% 3. Solve generalized eigenvalues problem

% Solve A*V = L*B*V (L = diag(lambda))
[V,L] = eigs(A,B,k,sigma);

% Extract eigenvalues as vector
eigval = diag(L);

% Extract physical eigenvectors
eigvec = V(1:ndof,:);

%% 4. Sort results

% Sort by proximity to sigma if sigma is a number
if isnumeric(sigma) && isscalar(sigma)
    % Compute squared distance
    dist = abs(eigval-sigma).^2;
    % Sort modes
    [~,id_sort] = sort(dist,'ascend');
    D = eigval(id_sort);
    V = eigvec(:,id_sort);
else
    D = eigval;
    V = eigvec;    
end

end