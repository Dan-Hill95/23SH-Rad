% Written by Dan J Hill (2021) - University of Surrey

% Inputs:    uu is the initial guess for a solution, 
%           p is all of the parameters in the problem,
%           mesh_params are the parameters of the mesh

% Outputs: F - Swift-Hohenberg system output for a given solution uu
% J - Jacobian of the function F.
function [F,J] = Equation_Patch(uu,p,mesh_params)

%Calling parameters
N  = mesh_params.N;         % Length of the mesh
LN = mesh_params.LN;        % Linear radial operator -(1+Laplacian)^2
LM = mesh_params.LM;        % D_{m} linear operator
u = uu(1:N);                % u - Amplitude of radial variable
v = uu(1+N:2*N);            % v - Amplitude of D_{m} variable

%% Pre-defining the size of functions, for speed
Fu = sparse(N,1);
Qu = sparse(N,1);
Cu = sparse(N,1);
Fv = sparse(N,1);
Qv = sparse(N,1);
Cv = sparse(N,1);

%% Explicit Forms
Lapu = LN - p(1)*speye(N);          % - (1 + d_rr + 1/r*d_r)^2 - mu
Lapv = LM - p(1)*speye(N);          % - (1 + d_rr + 1/r*d_r - (m/r)^2)^2 - mu
Qu(1:N) = u.^2 + 2.*v.^2;           % Quadratic nonlinearity for u equation
Cu(1:N) = (u.^2+6.*v.^2).*u;        % Cubic nonlinearity for u equation
Qv(1:N) = 2.*u.*v;                  % Quadratic nonlinearity for v equation
Cv(1:N) = 3.*(u.^2+v.^2).*v;        % Cubic nonlinearity for v equation


Fu = Lapu*u + p(2)*Qu + p(3)*Cu;    % Full u equation
Fv = Lapv*v + p(2)*Qv + p(3)*Cv;    % Full v equation
F = [Fu;Fv];                        % Output function F

if nargout > 1 %Calculating the Jacobian (only when called for)     
%% Pre-defining the size of functions, for speed

% Variations of u function with respect to u
Juu = sparse(N,N);
DQuu = sparse(N,N);
DCuu = sparse(N,N);
% Variations of u function with respect to v
Juv = sparse(N,N);
DQuv = sparse(N,N);
DCuv = sparse(N,N);
% Variations of v function with respect to u
Jvu = sparse(N,N);
DQvu = sparse(N,N);
DCvu = sparse(N,N);
% Variations of v function with respect to v
Jvv = sparse(N,N);
DQvv = sparse(N,N);
DCvv = sparse(N,N);

%% Explicit Forms 

% Variations of u function with respect to u
DQuu = sparse(1:N,[1:N], 2.*u ,N,N);
DCuu = sparse(1:N,[1:N], (u.^2+6.*v.^2) + 2.*u.^2 ,N,N);
Juu = Lapu + p(2)*DQuu + p(3)*DCuu;
% Variations of u function with respect to v
DQuv = sparse(1:N,[1:N], 4.*v ,N,N);
DCuv = sparse(1:N,[1:N], (12.*v).*u ,N,N);
Juv = p(2)*DQuv + p(3)*DCuv;
% Variations of v function with respect to u
DQvu = sparse(1:N,[1:N], 2.*v ,N,N);
DCvu = sparse(1:N,[1:N], 6.*(u).*v ,N,N);
Jvu = p(2)*DQvu + p(3)*DCvu;
% Variations of v function with respect to v
DQvv = sparse(1:N,[1:N], 2.*u ,N,N);
DCvv = sparse(1:N,[1:N], 3.*(u.^2+v.^2) + 6.*(v).*v ,N,N);
Jvv = Lapv + p(2)*DQvv + p(3)*DCvv;


J = [Juu, Juv; Jvu, Jvv];           %Output Jacobian J
end

end