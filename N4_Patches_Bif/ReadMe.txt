% Written by Dan Hill (2021) - University of Surrey, adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam; see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


Before running this code, make sure you have run the initialisation code

Init.m

in the previous folder, and have run the radial solver

Continuation_SH(p,SolnClass,Dir)

in the 'Radial_Solver' folder, for some choices of {p,SolnClass,Dir}. Then, there should be a data folder in 'Radial_Solver' called "[Type of initial guess]_Stab[m]_[Dir]" (where m is the last element of p).

% In order to run this code, you need to:

% 1. Call the function 'Cont_Patch_Bif'

% Example input: 
% Cont_Patch_Bif('pl');

% 2. Choose a data folder from 'Radial_Solver', and select the 'branch.mat' file.

% A bifurcation diagram will be shown, with linear D_{m} stability given by red/blue (unstable/stable) points. Each point has an 
% associated step number, for which every fifth is labelled. 

% 3. Choose a point that is either at the extreme of a stable/an unstable branch, and take note of its step number. Press any button to continue, % and a dialogue box will open. 

% 4. Choose the 'solution_[Step number].mat' file.


------------------------------------------------------------------
%% Input:   Cont_Patch_Bif(Dir)

%           SolnFolder              - This is chosen manually
%           File                    - This is chosen manually, should be close to a change in stability
%           Dir                     - Direction of parameter continuation;
%                                     must be either 'pl' (plus) or 'mn' (minus).


%% Purpose

% Solves a 4'th order Galerkin system for the stationary 2D 2-3 Swift-Hohenberg equation via finite-difference 
% methods and continues solutions in mu-parameter space. For U(r,theta) such that

% U(r,theta) = u[0](r) + 2*sum_{i=1}^{4} u[i](r)*cos(2*k*i*theta), 
% then each u[i](r) satisfies 

% 0= F[i](u[i]):= -(1+d^2_r + 1/r*d_r - (2*k*i/r)^2)^2 u[i] - mu*u[i] + f[i](U).

% Solving close to an initial guess for a standard radial profile, determined by SolnClass, perturbed by a destabilising D_{m} eigenmode,
% solutions are then continued in mu-parameter space.

%% Outputs

%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u[0]), L2Norm(u[1]), ..., L2Norm(u[4])]

% All data is stored in a folder named as "D[m]_Patch_[Solution Folder]_[Solution File]_[Dir]" 
