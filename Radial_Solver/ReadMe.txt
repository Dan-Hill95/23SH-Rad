% Written by Dan Hill (2021) - University of Surrey, adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam; see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


Before running this code, make sure you have run the initialisation code

Init.m

in the previous folder.

% In order to run this code, you only need to call the function 'Continuation_SH'

% Example input: 

% Continuation_SH([0.05, 1.6, -1, 6],'Sa','pl');

------------------------------------------------------------------
%% Input:   Continuation_SH(p,SolnClass,Dir)

%           p=[mu, nu, kappa, m]    - Initial parameters of system; 
%                                        mu is the bifurcation parameter, 
%                                        nu is the quadratic coefficient,
%                                     kappa is the cubic coefficient,
%                                         m is the dihedral order for computing stability,
%           SolnClass               - Type of radial solution for initial guess; 
%                                     must be either 'Sa' (Spot A), 'Sb' (Spot B), 'Ur' (Up-ring), or 'Dr' (Down-ring).
%           Dir                     - Direction of parameter continuation; must be either 'pl' (plus) or 'mn' (minus).

%% Purpose

% Solves the stationary axisymmetric 2-3 Swift-Hohenberg equation via finite-difference 
% methods and continues solutions in mu-parameter space. For u(r),

% 0= F(u):= -(1+d^2_r + 1/r*d_r)^2 u - mu*u + nu*u^2 + kappa*u^3.

% Solving close to an initial guess for a standard radial profile,
% determined by SolnClass, solutions are then continued by varying mu and solving again.
% At each step, linear stability is computed with respect to D_{m} perturbations.

%% Outputs

%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u)]

% All data is stored in a folder named as "[Type of initial guess]_Stab[m]_[Dir]" 
