% Adapted by Dan J Hill 2021
% Original code by David Lloyd (Surrey)
% Purpose: Setting up the radial grid and finite-difference matrices for
% solving the stationary axisymmetric Swift-Hohenberg equation

% m - is a dihedral lattice order, for computing linear stability with
% respect to D_{m} perturbations

%% Spatial grid

% Spatial coordinates: r direction
N = 1000; % Number of mesh points, must be even! 
L = 100; % domain truncation
h = L/(N-1); % space between each grid point
r = (0:N-1)'*h; % radial coordinate

%% Differentiation matrices

% compute radial Swift-Hohenberg operator differentiation matrix
e = ones(N,1);

% d_r
D = sparse(1:N-1,[2:N-1 N],ones(N-1,1)/2,N,N); 
D = (D - D')/(h);
D(1,2) = 0;
D(N,N-1) = 0;

% d_rr
D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],e,N,N);
D2 = (D2 + D2');

D2(1,2)=2;D2(N,N-1)=2;

D2 = D2/h^2;

% spy(D^2-D2);

r(1) = 1; 
R = sparse(1:N,[1:N],1./r,N,N); % 1/r*I (but well-defined at r=0)
r(1) = 0;
e = sparse(1:N,[1:N],e,N,N); % I

%% Radial Linear Operator

LN = D2 + R*D + e; % d_rr + d_r/r + I
LN(1,1) = 2*D2(1,1) + 1; LN(1,2) = 2*D2(1,2); % 2d_rr at r= 0
LN = -LN^2; % -(1+laplacian)^2

%% D_{m} Linear Operator for computing stability

LM = D2 + R*D + e - (m)^2*(R.^2).*e; % d_rr + d_r/r + I - (m/r)^2
LM(1,1) = (2-(m^2)/2)*D2(1,1) + 1; LM(1,2) = (2-(m^2)/2)*D2(1,2); % 2d_rr - (m^2/2) at r= 0
LM = -LM^2; % -(1+laplacian)^2

%% Saving everything in `mesh_params'

mesh_params.N  = N; 
mesh_params.L  = L; mesh_params.m = m; mesh_params.e = e;

mesh_params.r  = r; mesh_params.R = R;

mesh_params.D  = D; mesh_params.D2 = D2; mesh_params.LN = LN; mesh_params.LM = LM;


