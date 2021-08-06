%Written by Dan J Hill (2021) - University of Surrey

%Inputs:    uu is the initial guess for a solution, 
%           p is all of the parameters in the problem,
%           mesh_params are the parameters of the mesh

%Outputs:   F - Swift-Hohenberg equation output for a given solution uu
%           J - Jacobian of the function F.

function [F,J] = Equation_SH(uu,p,mesh_params)

%Calling parameters
N  = mesh_params.N;         % length of the mesh
LN = mesh_params.LN;        % Linear operator -(1+Laplacian)^2
u = uu(1:N);                % variable u

%Pre-defining the size of functions, for speed
F = sparse(N,1);
Q = sparse(N,1);
C = sparse(N,1);

Lap = LN - p(1)*speye(N);   % -(1+Laplacian)^{2} - mu
Q(1:N) = u.^2;              % Quadratic nonlinearity
C(1:N) = u.^3;              % Cubic nonlinearity


F = Lap*u + p(2)*Q + p(3)*C;    %Output function F

if nargout > 1 %Calculating the Jacobian (only when called for)
        
%Pre-defining the size of functions, for speed
J = sparse(N,N);
DQ = sparse(N,N);
DC = sparse(N,N);

DQ = sparse(1:N,[1:N],2.*u,N,N);        % Linearisation of the quadratic nonlinearity
DC = sparse(1:N,[1:N],3.*u.^2,N,N);     % Linearisation of the cubic nonlinearity
J = Lap + p(2)*DQ + p(3)*DC;            % Output Jacobian J
end
end