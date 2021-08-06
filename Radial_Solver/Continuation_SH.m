% Written by Dan Hill (2021) - University of Surrey, 
% adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam;
% see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.

% Example input: 
% Continuation_SH([0.05, 1.6, -1, 6],'Sa','pl');

%% Inputs
%           p=[mu, nu, kappa, m]    - Initial parameters of system; 
%                                     mu is the bifurcation parameter, 
%                                     nu is the quadratic coefficient,
%                                     kappa is the cubic coefficient,
%                                     m is the dihedral order for computing stability,
%           SolnClass               - Type of radial solution for initial guess; 
%                                     must be either 'Sa' (Spot A), 'Sb' (Spot B), 'Ur' (Up-ring), or 'Dr' (Down-ring).
%           Dir                     - Direction of parameter continuation;
%                                     must be either 'pl' (plus) or 'mn' (minus).

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

function branch = Continuation_SH(p,SolnClass,Dir)

close all, clc;

%% Initial Set-up
%Define mesh parameters
m=p(4);             % Need to assign m to p(4) prior to running SetupDiffMats_SH
if mod(m,1)==1
    error('m must be a natural number.');
end
SetupDiffMats_SH;   % Defines mesh parameters as 'mesh_params'

% initial guess
p0=p;
[u0,SolnName]=InitialGuess_SH(p0,SolnClass,mesh_params);

% Plotting initial data
hfig1 = figure;
 t = 0:0.01:2*pi;
   t = t';
  [R,T]=meshgrid(r,t);
  [U0,T0]=meshgrid(u0,t);
  
  % contour plot
  subplot(1,2,1)  
  z=surf(R.*cos(T), R.*sin(T),U0);
  view(45,25);
  pbaspect([1 1 1]);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
  axis([-30 30 -30 30 -2*max(u0) 3*max(u0)]);
  
  % radial profile
  subplot(1,2,2)
  plot(r,u0,'.-');
  pbaspect([1 1 1]);
  axis([0 L -1.2*max(-u0) 1.2*max(u0)]);
pause;    

%% Converging initial guess
    myproblemHandle = @(u) Equation_SH(u,p0,mesh_params);
    
options = optimset('Jacobian','on','Display','iter','MaxIter',30,'TolFun',1e-7,'DerivativeCheck','off');

%Solving the radial Swift-Hohenberg equation close to the initial guess

[u_out,fval,exitflag,output,jacobian] = fsolve(myproblemHandle,u0,options);

%Plotting the converged solution

hfig1; % rewriting previous plot, for space

u2 = u_out(1:N);
[U,T0]=meshgrid(u2,t);

% contour plot
subplot(1,2,1)  
z=surf(R.*cos(T), R.*sin(T),U);
  view(45,25);
  pbaspect([1 1 1]);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
  axis([-30 30 -30 30 -2*max(u2) 3*max(u2)]);
% radial profile 
  subplot(1,2,2)
  plot(r,u2(1:N),'.-');
  pbaspect([1 1 1]);
  axis([0 L -1.2*max(-u2) 1.2*max(u2)]);
  pause;
  
  close(hfig1); % Closing plots prior to continuation, for space
  
%% Define handle to right-hand side and time output function
prob = @(u,p) Equation_SH(u,p,mesh_params);
plotSol  = @(u,p,parent) PlotSurface_SH(u,parent,mesh_params);
solMeas  = @(step,u,p) SolutionMeasures_SH(u,mesh_params);
compSpec = @(u,p) ComputeEigenvalues_SH(u,p,mesh_params);
plotSpec = @(d,p,parentHandle) PlotSpectrum_SH(d,p,parentHandle);

%% Assign problem 
stepperPars.iContPar      = 1;
% Choosing which direction to perform secant continuation
if Dir == 'pl'
stepperPars.s0            = 0.01;
elseif Dir== 'mn' 
stepperPars.s0            = -0.01;
else
    error('Final argument must be "pl" or "mn".');
end
stepperPars.sMin          = 1e-8;
stepperPars.sMax          = .1;
stepperPars.pMin          = 0;
stepperPars.pMax          = 2;
stepperPars.maxSteps      = 1000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'TolFun',1e-4,...
                                     'Jacobian','on',...
                                     'MaxIter',30);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = [SolnName '_Stab' num2str(m) '_' Dir ];
stepperPars.PlotSolution  = plotSol;
stepperPars.BranchVariables = solMeas;
stepperPars.ComputeEigenvalues = compSpec;
stepperPars.PlotSpectrum = plotSpec;      
stepperPars.PlotBranchVariableId = 1;


branch = SecantContinuation(prob,u_out,p0,stepperPars,'Branch');
%
end