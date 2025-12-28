function [N,dN,d2N,d3N] = getShapeFun(x)
%GETSHAPEFUN Evaluates Cubic Hermite shape functions and derivatives.
%
%   Computes the shape functions for a 2-node, 4-DOF reference element
%   defined on xi in [-1, 1].
%   DOFs: [Node1_Value, Node1_Slope, Node2_Value, Node2_Slope]
%
%   Inputs:
%       x   : Coordinate in element domain [-1, 1].
%
%   Outputs:
%       N   : [4 x 1] Shape function values.
%       dN  : [4 x 1] First derivatives (dN/dxi).
%       d2N : [4 x 1] Second derivatives (d^2N/dxi^2).
%       d3N : [4 x 1] Third derivatives (d^3N/dxi^3).

% Shape functions (N) (cubic) 
N = [
    0.25*( 2 - 3*x        + x.^3);
    0.25*( 1 -   x - x.^2 + x.^3);
    0.25*( 2 + 3*x        - x.^3);
    0.25*(-1 -   x + x.^2 + x.^3);
];

% First derivatives (dN/dx) (quadratic)
dN = [
    0.25*(-3 +     + 3*x.^2);
    0.25*(-1 - 2*x + 3*x.^2);
    0.25*( 3 +     - 3*x.^2);
    0.25*(-1 + 2*x + 3*x.^2);
];

% Second derivatives (d^2N/dx^2) (linear)
d2N = [
    0.25*(     6*x);
    0.25*(-2 + 6*x);
    0.25*(   - 6*x);
    0.25*( 2 + 6*x);
];

% Third derivatives (d^3N/dx^3) (constant)
if nargout > 3
    d3N = [
        1.5;
        1.5;
       -1.5;
        1.5;
    ].*ones(size(x));
end

end