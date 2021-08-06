% Written by Dan Hill (2021) - University of Surrey, adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam; see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


Before running this code, make sure you have run the initialisation code

Init.m

in the previous folder, and have run the radial solver

Continuation_SH(p,SolnClass,Dir)

in the 'Radial_Solver' folder, for some choices of {p,SolnClass,Dir}. Then, there should be a data folder in 'Radial_Solver' called "[Type of initial guess]_Stab[m]_[Dir]" (where m is the last element of p).

% In order to run this code, you need to:

% 1. Call the function 'Continuation_Bif'

% Example input: 
% Continuation_Bif(2,'mn');

% 2. Choose a data folder from 'Radial_Solver', and select the 'branch.mat' file.

% A bifurcation diagram will be shown, with linear D_{m} stability given by red/blue (unstable/stable) points. Each point has an 
% associated step number, for which every fifth is labelled. 

% 3. Choose a point that is either at the extreme of a stable/an unstable branch, and take note of its step number. Press any button to continue, % and a dialogue box will open. 

% 4. Choose the 'solution_[Step number].mat' file.


------------------------------------------------------------------
%% Input:   Continuation_Bif(ContVar,Dir)

%           ContVar                 - Denotes which parameter value in p will be varied during continuation
%           SolnFolder              - This is chosen manually
%           File                    - This is chosen manually, should be close to a change in stability
%           Dir                     - Direction of parameter continuation; must be either 'pl' (plus) or 'mn' (minus).

%% Purpose

% Solves the linear stability problem for the stationary axisymmetric 2-3 Swift-Hohenberg equation close to a bifurcation point for some
% localised radial solution, and then continues the bifurcation point through parameter space. For u(r) that solves

% 0 = F(u) := -(1+d^2_r + 1/r*d_r)^2 u - mu*u + nu*u^2 + kappa*u^3,

% we look for a destabilising D_{m} eigenmode v(r) such that,

% 0 = J[u]*v := [-(1+d^2_r + 1/r*d_r - (m/r)^2)^2 - mu + 2*nu*u + 3*kappa*u^2]*v,
% 0 = |v|^2-1.

% An initial solution [File] is chosen from the data folder [SolnFolder] close to a change in D_{m} stability. 
% Then, the bifurcation point is tracked through parameter space as the parameters mu and p(ContVar) are varied.

% Here, we recall that  p = [mu, nu, kappa, m], and so we only really want ContVar = 2 or 3.

%% Outputs

%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u), L2Norm(v), p(ContVar)]

% All data is stored in a folder named as "Bif_[Solution Folder]_[Solution File]_[Dir]" 