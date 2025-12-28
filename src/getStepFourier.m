function [H_,k_] = getStepFourier(a,b,shape,N,sym)
%GETSTEPFOURIER Computes Fourier coefficients for a 2D periodic step function.
%
%   [H_, k_] = GETSTEPFOURIER(a, b, shape, N, sym) calculates the N most relevant
%   Fourier coefficients (H) and associated wavenumbers (k) for a PSub 
%   (periodic substrate) shape defined by a step function.
%
%   The domain is defined as:
%         ___________
%        /          /
%       /  ___     /     
%   a2 /  |   |   /          
%     /   |___|  /             
%    /    <-b-> /               
%   /__________/               
%       a1
%
%   Inputs:
%       a     : [a1, a2] Periodicity lengths in each direction.
%       b     : Width of the step function (PSub) at the center.
%       shape : 'none', 'fullspan', 'square', or 'hexagon'.
%       N     : Number of modes to retain (sorted by magnitude).
%       sym   : (Boolean) If true, exploits spanwise symmetry (removes k_y < 0).
%
%   Outputs:
%       H_    : [N x 1] Fourier coefficients.
%       k_    : [N x 2] Fourier wavenumbers [alpha, beta].

%% 1. Setup and pre-calculations

% Check inputs
if isscalar(a), a = [a,a]; end

% Define PSub dimensions
if strcmpi(shape,'fullspan')  || strcmpi(shape,'none')
    b = [b,a(2)/2]; % The PSub covers the entire span
else
    b = [b,b];
end

%% 2. Generate candidate modes

% Get a grid large enough to find the best N modes
search_range = max(50, ceil(sqrt(N)*3));
n_ = getFourierIndices(search_range,'full');

% Map integer indices to physical wavenumbers
[k_,a_] = getFourierWavenumbers(n_,a,shape);

%% 3. Filter modes

% If the 'sym' option is activated, remove redundant modes
if sym
    k_(k_(:,2)<0,:) = [];
end

% For 1D-like shapes remove spanwise variations
if strcmpi(shape,'fullspan')  || strcmpi(shape,'none')
    k_(abs(k_(:,2))>0,:) = [];
end

%% 4. Compute Fourier coefficients

% Normalization constant
c = 1/norm(cross([a_(1,:),0],[a_(2,:),0]));

% Streamwise step
kx = k_(:,1);
idx = abs(kx)>0;
A = 2*b(1)*ones(size(kx)); % Limit as kx -> 0
A(idx) = sin(kx(idx)*b(1))./(kx(idx)/2);

% Spanwise step
kz = k_(:,2);
idz = abs(kz)>0;
B = 2*b(2)*ones(size(kz)); % Limit as kz -> 0
B(idz) = sin(kz(idz)*b(2))./(kz(idz)/2);

% Amplitude coefficients
H_ = c*(A.*B);

% Sort modes according to relevance 
[~,ind] = sort(abs(H_),'descend');

% Retain N most relevant modes
lim = min([N,length(ind)]);
H_ = H_(ind(1:lim));
k_ = k_(ind(1:lim),:);

end

%--------------------------------------------------------------------------
% LOCAL HELPER FUNCTIONS
%--------------------------------------------------------------------------

function [n_] = getFourierIndices(N,opt)
%GETFOURIERINDICES Generates a grid of integer pairs (n, m).

% Check inputs
if isscalar(N)
    Nx = N; Nz = N;
else
    Nx = N(1); Nz = N(2);
end

% Create initial grid (upper half plane)
i_ = [0,1:Nx,repmat(-Nx:Nx,1,Nz)]';
j_ = [zeros(1,Nx+1),repelem(1:Nz,1,2*Nx+1)]'; 
n_ = [i_(:),j_(:)];

% Extend indices to full plane (if requested)
if exist('opt','var') && strcmpi(opt,'full')
    % Mirror indices to get negative j components
    n_ = [n_(1,:);reshape([-n_(2:end,:)';n_(2:end,:)'],2,[])'];
end

end

function [k_,a_] = getFourierWavenumbers(n_,a,shape)
%GETFOURIERWAVENUMBERS Maps integer indices to physical wavenumbers.

% Determine rotation angles for periodicity vectors
if isa(shape,'char') % If shape is given
    switch shape
        case {'square','fullspan','none'}
            ca = 1; sa = 0;
            cb = 1; sb = 0;
        case 'hexagon' % Hexagonal lattice (60 deg)
            ca = sqrt(3)/2; sa = 1/2;
            cb = 1;         sb = 0;
    end
else % Custom angle input
    ca = cos(shape(1)); sa = sin(shape(1));
    cb = cos(shape(2)); sb = sin(shape(2));
end

% Determine the periodicity (lattice) vectors
a_ = [a(1)*ca, a(1)*sa;   % vector a1
     -a(2)*sb, a(2)*cb];  % vector a2

% Compute associated wavenumbers
det = ca*cb+sa*sb;
k_ = zeros(size(n_));
k_(:,1) = (2*pi/det)*(n_(:,1)*cb/a(1)-n_(:,2)*sa/a(2));
k_(:,2) = (2*pi/det)*(n_(:,1)*sb/a(1)+n_(:,2)*ca/a(2));

end