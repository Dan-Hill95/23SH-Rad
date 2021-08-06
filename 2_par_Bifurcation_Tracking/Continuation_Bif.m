% Written by Dan Hill (2021) - University of Surrey, 
% adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam;
% see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


%% Inputs
%           ContVar                 - Denotes which parameter value in p will be varied during continuation
%           SolnFolder              - This is chosen manually
%           File                    - This is chosen manually, should be close to a change in stability
%           Dir                     - Direction of parameter continuation;
%                                     must be either 'pl' (plus) or 'mn' (minus).

%% Purpose
% Solves the linear stability problem for the stationary axisymmetric 2-3
% Swift-Hohenberg equation close to a bifurcation point for some localised radial solution, and then continues
% the bifurcation point through parameter space. For u(r) that solves

% 0 = F(u) := -(1+d^2_r + 1/r*d_r)^2 u - mu*u + nu*u^2 + kappa*u^3,

% we look for a destabilising D_{m} eigenmode v(r) such that,

% 0 = J[u]*v := [-(1+d^2_r + 1/r*d_r - (m/r)^2)^2 - mu + 2*nu*u + 3*kappa*u^2]*v,
% 0 = |v|^2-1.

% An initial solution [File] is chosen from the data folder [SolnFolder]
% close to a change in D_{m} stability. Then, the bifurcation point is
% tracked through parameter space as the parameters mu and p(ContVar) are varied.

% Here, we recall that  p = [mu, nu, kappa, m], and so we only really want
% ContVar = 2 or 3.

%% Outputs
%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u), L2Norm(v), p(ContVar)]

% All data is stored in a folder named as "Bif_[Solution Folder]_[Solution File]_[Dir]" 

function branch = Continuation_Bif(ContVar,Dir)

close all, clc;

%% Initial Set-up

K = pwd;
cd ..
cd Radial_Solver
% Select which radial solution to investigate
[fileName0,pathName] = uigetfile('*.mat','Select the branch file');
branchFile = [pathName fileName0];
cd(pathName)
ExploreBifurcationDiagram(branchFile,5);
% Choose a point close to a bifurcation
pause;
[fileName,pathName] = uigetfile('*.mat','Select the solution file');
F = [pathName fileName];
% Defining the folder name
parts=strsplit(pathName,'\');
solnFolder=char(parts(end-1));
% Defining the file number
parts1=strsplit(fileName,'.');
parts2=strsplit(char(parts1(1)),'_');
File=str2num(char(parts2(2)));
sol = load(F);
close all
% Returning to original directory
cd(K)
p = sol.p;                                                  % p is determined from the initial data
m=p(4);                                                     % m must be labelled before running SetupDiffMats_SH
if mod(m,1)==1
    error('m must be a natural number.');
end
uu0= sol.u;
SetupDiffMats_SH;                                           % Defines mesh parameters as 'mesh_params'
mesh_params.ContVar=ContVar;

my_rhs = @(u) Equation_SH(u,p,mesh_params);                 % Checks that our data is a solution of the SH Equation

options = optimset('Jacobian','on','Display','iter','MaxIter',1000,'DerivativeCheck','off','Algorithm','levenberg-marquardt');

[u_out,fval,exitflag,output,jacobian] = fsolve(my_rhs,uu0,options);

% Plotting initial solution
r=mesh_params.r;
N=mesh_params.N;
w_out = u_out(1:N);
hfig1=figure;
plot(r,w_out);
pause;
close(hfig1);

%% Linear Stability Problem

%Compute the eigenmode V for the 0 eigenvalue
[V,~] = eigs(jacobian,1,0.01);
myproblemHandle = @(u,p)  Equation_SH(u,p,mesh_params);
my_rhs2 = @(u) bif_cont_SH(u,p,mesh_params,myproblemHandle); % Set up the 2-parameter system for tracking bifurcations

% Define an initial guess of the original solution, the destabilising
% eigenmode, and the parameter that will be varied
uu1 = [u_out; V(:,1); p(ContVar)];
options = optimset('Jacobian','on','Display','iter','MaxIter',20,'TolFun',1e-4,'DerivativeCheck','off');

[u_out2,fval,exitflag,output,jacobian] = fsolve(my_rhs2,uu1,options);


%% Define handle to right-hand side and time output function
problemHandle            = @(u,p)  bif_cont_SH(u,p,mesh_params,myproblemHandle);
plotSolutionHandle       = @(u,p,parentHandle) PlotSolution_fold(u,parentHandle,mesh_params);
branchVariablesHandle    = @(step,u,p) SolutionMeasures_fold(u,mesh_params);
computeEigenvaluesHandle = [];%@(u,p) ComputeEigenvalues(u,p,t);
plotSpectrumHandle       = []; %@(d,p,parentHandle) PlotSpectrum(d,p,parentHandle);
% 
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
stepperPars.pMin          = -2;
stepperPars.pMax          = 2;
stepperPars.maxSteps      = 100;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'TolFun',1e-4,...
                                     'Jacobian','on',...
                                     'MaxIter',20);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = ['Bif_' solnFolder '_' num2str(File) '_' Dir];
stepperPars.PlotSolution  = plotSolutionHandle;
stepperPars.BranchVariables = branchVariablesHandle;
stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = computeEigenvaluesHandle;
stepperPars.PlotSpectrum = plotSpectrumHandle;      
stepperPars.PlotBranchVariableId = 1;

branch = SecantContinuation(problemHandle,u_out2,p,stepperPars,'FoldCont');


