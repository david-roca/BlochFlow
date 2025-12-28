%% Example usage: PSub-channel flow Bloch analysis
clear; clc;
addpath(genpath('src')); % Automatically finds the functions

%% 1. Define parameters

% Fluid properties
fluid.Re = 7500;         % Reynolds number
fluid.delta = 4.38e-4;   % Reference length, half channel height (m) 
fluid.rho = 1000;        % Density (kg/m3)
fluid.mu = 1e-3;         % Viscosity (Pa·s)
fluid.shape = 'hexagon'; % Options: 'none', 'fullspan', 'square', 'hexagon'
% Calculate reference velocity
fluid.Uc = fluid.Re*fluid.mu/fluid.rho/fluid.delta;
% Calculate reference TS mode wavelength
fluid.lambda_TS = 2*pi*fluid.delta;
% Calculate fluid unit cell size
fluid.a = fluid.lambda_TS/2;

% Solid/PSub properties
solid.Ncell = 5;         % Number of unit cells
solid.lrod = 0.01;       % Unit cell thickness (m)
solid.E = 3e9;           % Young's modulus (Pa)
solid.rho = 1200;        % Solid density (kg/m3)
solid.q = [0,6e-8];      % Damping factor
solid.mres = 10;         % Resonators mass ratio
solid.fres = 2000;       % Resonators frequency (Hz)
solid.e = 0.5;           % PSub size ratio
% Calculate PSub cross-section area
solid.b = solid.e*fluid.a/2;
if strcmpi(fluid.shape,'fullspan') || strcmpi(fluid.shape,'none')
    solid.A = 2*solid.b*fluid.a;
else
    solid.A = (2*solid.b)^2;
end
% Calculate resonators mass
solid.mr_ = solid.mres*solid.rho*solid.lrod*solid.A;
% Calculate resonators stiffness
solid.kr_ = (2*pi*solid.fres)^2*solid.mr_;

% Numerical settings
num.Ny = 100;            % Number of elements in the vertical direction
num.Nx = 13;             % Number of Fourier modes
num.Nm = 5;              % Number of Bloch modes
num.k_target = 1;        % Guessed wavenumber (non-dimensional)

%% 2. Run solver

% Compute solid matrices
[K,M,C] = computeSolidMatrices(solid.lrod,solid.A,solid.rho,solid.E, ...
                               solid.q,solid.mr_,solid.kr_,solid.Ncell);

% Test frequency
freq = 1555; % Hz

% Compute admittance
if strcmpi(fluid.shape,'none')
    YPSub = 0;
else
    YPSub = computeAdmittance(2*pi*freq,K,M,C,solid.A);
end

% Non-dimensional frequency
w = 2*pi*freq/(fluid.Uc/fluid.delta);

% Non-dimensional admittance
YPSub_nd = fluid.rho*fluid.Uc*YPSub;

% Non-dimensional unit cell size
a_nd = fluid.a/fluid.delta;
b_nd = solid.b/fluid.delta;

disp(['Starting computation at Re = ',num2str(fluid.Re), ...
      ', Y = ',num2str(YPSub_nd),...
      ', and f = ',num2str(freq),' Hz...']);

% Solve Bloch wave problem in the fluid
[k,v_] = solveBlochFluid(w,num.k_target,fluid.Re,YPSub_nd,a_nd,b_nd,...
                         fluid.shape,num.Nx,num.Ny,num.Nm,'cosine');

%% 3. Display results
disp(['Computation complete! Resulting wavenumber (k): ', num2str(k(1))]);