% Written by Dan Hill (2021) - University of Surrey, 
% adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam;
% see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.


%% Inputs
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

% Solving close to an initial guess for a standard radial profile,
% determined by SolnClass, perturbed by a destabilising D_{m} eigenmode,
% solutions are then continued in mu-parameter space.

%% Outputs
%           branch- [Step, Stability, mu, EucNorm(u), L2Norm(u[0]), L2Norm(u[1]), ..., L2Norm(u[4])]

% All data is stored in a folder named as "D[m]_Patch_[Solution Folder]_[Solution File]_[Dir]" 

function branch = Cont_Patch_Bif(Dir)

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
p=sol.p;
m = p(4);  
if mod(m,2)>0
    error('m must be an even natural number.');
end
SetupDiffMats_4Patch;
ContVar=2;
mesh_params.ContVar=ContVar;

% initial data
N = mesh_params.N;
r = mesh_params.r;
uu0=sol.u;

%% Initial Perturbation
uu1=uu0;
uv1=0.2*norm(uu1,2)*((-1)^m)*besselj(m,r).*exp(-sqrt(p(1))*r/2);
hfig1=figure;
t = 0:0.01:2*pi;
   t = t';
  u2 = [uu1;uv1];
  [R,T]=meshgrid(r,t);
  [U,T0]=meshgrid(u2,t);
  UU(:,:,1)=U(:,1:N)/2;
  for i=1
      UU(:,:,i+1)=U(:,1+i*N:(i+1)*N);
      UU(:,:,i+1)=UU(:,:,i+1).*cos(m.*T);
  end
  z=surf(R.*cos(T), R.*sin(T),sum(UU,3));
  view(45,25);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
  axis([-30 30 -30 30 -1.2*max(-uu1) 1.2*max(uu1)]);
  pause;
close(hfig1);
u1=[uu1;uv1;zeros(3*N,1)]; 
%% Define handle to right-hand side and time output function
prob     = @(u,p) Swift_Patch_N4(u,p,mesh_params);
plotSol  = @(u,p,parent) PlotSurface_Patch(u,p,parent,mesh_params);
solMeas  = @(step,u,p) SolutionMeasures_Patch(step,u,p,mesh_params);
compSpec = @(u,p) ComputeEigenvalues_Patch(u,p,prob);
plotSpec = @(d,p,parentHandle) PlotSpectrum_SH(d,p,parentHandle);

%% Assign problem 
stepperPars.iContPar      = 1;
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
stepperPars.maxSteps      = 200;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'TolFun',1e-4,...
                                     'Jacobian','on',...
                                     'MaxIter',30);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = ['D',num2str(m),'_Patch_' solnFolder '_' num2str(File) '_' Dir];
stepperPars.PlotSolution  = plotSol;
stepperPars.BranchVariables = solMeas;
% stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = compSpec;
stepperPars.PlotSpectrum = plotSpec;      
stepperPars.PlotBranchVariableId = 1;


branch = SecantContinuation(prob,u1,p,stepperPars,'Patch_Cont1');
%
end