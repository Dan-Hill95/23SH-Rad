% Written by Dan Hill (2021) - University of Surrey, adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam; see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


Before running this code, make sure you have run the initialisation code

Init.m

in the previous folder, and have run the radial solver

Continuation_SH(p,SolnClass,Dir)

in the 'Radial_Solver' folder, for some choices of {p,SolnClass,Dir}. Then, there should be a data folder in 'Radial_Solver' called "[Type of initial guess]_Stab[m]_[Dir]" (where m is the last element of p).

% In order to run this code, you need to:

% 1. Call the function 'Continuation_Patch'

% Example input: 
% Continuation_Patch(0.1,'pl');

% 2. Choose a data folder from 'Radial_Solver', and select the 'branch.mat' file.

% A bifurcation diagram will be shown, with linear D_{m} stability given by red/blue (unstable/stable) points. Each point has an 
% associated step number, for which every fifth is labelled. 

% 3. Choose a point that is either at the extreme of a stable/an unstable branch, and take note of its step number. Press any button to continue, % and a dialogue box will open. 

% 4. Choose the 'solution_[Step number].mat' file.


------------------------------------------------------------------
%% Input:   Continuation_Patch(Pert,Dir)
%           Pert                    - Denotes the magnitude and sign of the
%                                     D_{m} perturbation
%           SolnFolder              - This is chosen manually
%           File                    - This is chosen manually, should be close to a change in stability
%           Dir                     - Direction of parameter continuation;
%                                     must be either 'pl' (plus) or 'mn' (minus).


%% Purpose
% Solves a coupled Galerkin system for the stationary 2D 2-3 Swift-Hohenberg equation via finite-difference 
% methods and continues solutions in mu-parameter space. For u(r) and v(r) such that

% U(r,theta) = u(r) + 2*v(r)*cos(m*theta), then u(r) and v(r) satisfy 

% 0= F1(u,v):= -(1+d^2_r + 1/r*d_r)^2 u - mu*u + nu*[u^2 + 2*v^2] + kappa*[u^2 + 6*v^2]*u,
% 0= F2(u,v):= -(1+d^2_r + 1/r*d_r - (m/r)^2)^2 v - mu*v + 2*nu*u*v + 3*kappa*[u^2 + v^2]*v.

% Solving close to an initial guess for a standard radial profile,
% determined by SolnClass, perturbed by a destabilising D_{m} eigenmode,
% solutions are then continued in mu-parameter space.

%% Outputs
%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u), L2Norm(v)]

% All data is stored in a folder named as "Patch_[Solution Folder]_[Solution File]_[Dir]" 
