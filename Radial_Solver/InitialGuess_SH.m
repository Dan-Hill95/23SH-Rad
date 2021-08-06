% Written by Dan J Hill 2021
%% Inputs
%           p               - Parameters of the system
%           SolnClass       - Type of radial pattern for the initial guess
%           mesh_params     - Parameters of the mesh

%% Purpose
% To solve the stationary axisymmetric 2-3 Swift-Hohenberg equation close
% to some standard localised pattern (Spots A/B, or Up/Down Rings)

%% Outputs
%           uout            - Converged localised radial solution
%           SolnName        - Type of radial pattern, SolnClass converted into a full name

function [uout,SolnName] = InitialGuess_SH(p,SolnClass,mesh_params)
% Set-up r-coordainte and fsolve options
r=mesh_params.r;
options=optimset('Display','iter','Jacobian','on','MaxIter',100);     % Option to display output and use Jacobian

%% Initial profiles for standard localised radial patterns

if SolnClass == 'Sa'            %Spot A
u=sqrt(p(1)).*(sqrt(3)/p(2))*besselj(0,r).*exp(-sqrt(p(1)).*r./2);
SolnName= 'SpotA';
elseif SolnClass == 'Sb'        %Spot B
u=-0.5*besselj(0,r).*exp(-sqrt(p(1)).*r./2);
SolnName= 'SpotB';
elseif SolnClass == 'Ur'        %Up Ring
u=0.1.*r.*besselj(1,r).*exp(-sqrt(p(1)).*r./2).*(1-exp(-sqrt(p(1)).*r./2));
SolnName= 'Up_Ring';
elseif SolnClass == 'Dr'        %Down Ring
% Note: I've had more success finding down ring solutions by first finding
% up rings and then taking the mirror image.
uu0=0.1.*r.*besselj(1,r).*exp(-sqrt(p(1)).*r./2).*(1-exp(-sqrt(p(1)).*r./2)); %Up ring guess
[uu,fval,exitflag,output,jac_fsolve] = fsolve(@(u) Equation_SH(u,p,mesh_params),uu0,options);
 t = 0:0.01:2*pi;
   t = t';
  [R,T]=meshgrid(r,t);
  [U0,T0]=meshgrid(uu,t);
% Contour plot
subplot(1,2,1)         
z=surf(R.*cos(T), R.*sin(T),U0);
  view(0,90);
  pbaspect([1 1 1]);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
  axis([-50*sqrt(2) 50*sqrt(2) -50*sqrt(2) 50*sqrt(2) -2*max(uu) 3*max(uu)]);
% Radial Profile
  subplot(1,2,2)
  plot(r,uu,'.-');
  pbaspect([1 1 1]);
  axis([0 100 -1.2*max(-uu) 1.2*max(uu)]);
pause;    
close all
% Define down ring as the mirror image of an up ring
u=-uu;
SolnName= 'Down_Ring';
else
    error('The second argument must be one of the following solution classes: "Sa" (Spot A), "Sb" (Spot B), "Ur" (Up Ring), or "Dr" (Down Ring).');
end
u0=u;

% Solve the Swift-Hohenberg equation close to the choen standard localised radial pattern
[uout,fval,exitflag,output,jac_fsolve] = fsolve(@(u) Equation_SH(u,p,mesh_params),u0,options);
end