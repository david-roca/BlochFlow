function [U,dU,d2U] = getBaseFlow(y)
%GETBASEFLOW Evaluates the base flow velocity profile and its derivatives.
%
%   Computes the non-dimensional Plane Poiseuille flow profile for a channel 
%   defined between y=0 (bottom wall) and y=2 (top wall).
%   
%   Profile: U(y) = 1 - (1-y)^2  (Parabolic)
%   
%   Inputs:
%       y   : [N x M] Array of wall-normal coordinates (0 <= y <= 2).
%
%   Outputs:
%       U   : Base flow velocity at y.
%       dU  : First derivative (dU/dy) (Shear rate).
%       d2U : Second derivative (d^2U/dy^2) (Curvature).

% Base profile: Parabola with max value 1 at y=1, and 0 at walls
U = 1-(1-y).^2;

% First derivative: 2*(1-y)
% Positive at bottom wall (y<1), zero at center, negative at top.
dU = 2*(1-y);

% Second derivative: Constant -2
d2U = -2*ones(size(y));
    
end