function [uout,r] = InitialGuess_Patch(p,mesh_params)
% clear all
close all

m=p(4);

N = mesh_params.N;
m = mesh_params.m;
r = mesh_params.r;
mesh_params1=mesh_params;
LM = mesh_params.LM;
mesh_params1.LM = LM(1:N,1:N);

% initial data
C = 2*cos(m*pi/3);

ba=1/C;

Ba=sqrt((C-1)/(C^3));

u00 =0.5*(sqrt(3)/p(2))*(ba).*besselj(0,r).*exp(-sqrt(p(1)).*r./2); 
v00= zeros(N,1);
v00(1:N) = 0.5*(sqrt(3)/p(2))*Ba*((-1)^(m))*besselj(m,r).*exp(-sqrt(p(1)).*r./2);    
u0=[u00; v00];


options=optimset('Display','iter','Jacobian','on','MaxIter',100);     % Option to display output and use Jacobian

[uout,fval,exitflag,output,jac_fsolve] = fsolve(@(u) Equation_Patch(u,p,mesh_params1),u0,options);              % Call fsolve

end