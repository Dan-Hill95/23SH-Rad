% Written by Dan Hill (2021) - University of Surrey, adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam; see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


Before running this code, make sure you have run the initialisation code

Init.m

in the previous folder.

% In order to run this code, you only need to call the function 'Cont_Patch_N1'

% Example input: 

% Cont_Patch_N1([0.02, 1.6, -1, 2],'pl');

------------------------------------------------------------------
%% Input:   Cont_Patch_N1(p,'pl')

%           p=[mu, nu, kappa, m]    - Initial parameters of system; 
%                                        mu is the bifurcation parameter, 
%                                        nu is the quadratic coefficient,
%                                     kappa is the cubic coefficient,
%                                         m is the dihedral lattice of the solution (must be even),
%           Dir                     - Direction of parameter continuation; must be either 'pl' (plus) or 'mn' (minus).


%% Purpose

% Solves a 4'th order Galerkin system for the stationary 2D 2-3 Swift-Hohenberg equation via finite-difference 
% methods and continues solutions in mu-parameter space. For U(r,theta) such that

% U(r,theta) = u[0](r) + 2*sum_{i=1}^{4} u[i](r)*cos(m*i*theta), then each u[i](r) satisfies 

% 0= F[i](u[i]):= -(1+d^2_r + 1/r*d_r - (m*i/r)^2)^2 u[i] - mu*u[i] + f[i](U).

% Solving for a small truncated patch obtained analytically, solutions are then continued in mu-parameter space.

%% Outputs

%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u[0]), L2Norm(u[1]), ..., L2Norm(u[4])]

% All data is stored in a folder named as "D[m]_Patch_[Solution Folder]_[Solution File]_[Dir]" 