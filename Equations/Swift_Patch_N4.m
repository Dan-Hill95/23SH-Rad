%Written by Dan Hill (2020) - University of Surrey

%Inputs: uu is all of our solutions, p is all of our parameters in the problem,
%mesh_params are the parameters of the mesh

%Coupled system of 4'th-order Galerkin truncation of the 2D
%Swift-Hohenberg equation in polar coordinates. For u(r,theta),

%0= F(u):= -(1+d^2_r + 1/r*d_r + 1/r^2*d^2_theta)^2 u - mu*u + a1*u^2 + a4*u^3.

% Solutions take the form u(r,theta) = u0(r)/2 + sum_{i=1}^{4} ui(r)*cos(2*k*i*theta)

%By projecting onto each Fourier mode, we arrrive at 5 equations to
%solve, found below.

%Outputs: F - Coupled Swift-Hohenberg system output for a given solution
%uu, J - Jacobian of the function F.
function [F,J] = Swift_Patch_N4(uu,p,mesh_params)

%Calling parameters
N = mesh_params.N;%length of the mesh
LN = mesh_params.LN;%-(1+d^2_r + 1/r*d_r)^2 finite difference matrix
LM = mesh_params.LM;%-(1+d^2_r + 1/r*d_r - (m*i)^2/r^2)^2 finite difference matrix for each u[i](r)


u = uu(1:5*N);%uu holds information on the size of the perturbations. We don't want that information when using u, hence we omit those components of uu

u0 = u(1:N);
u1 = u(1+N:2*N);
u2 = u(1+2*N:3*N);
u3 = u(1+3*N:4*N);
u4 = u(1+4*N:5*N);

%renaming parameters coming from p
mu = p(1);

%Defining the linear operator L:=-(1+d^2_r + 1/r*d_r - (2*k*m)^2/r^2)^2 - mu,
%for each m = 1, 2,..., x.
Lap(1:N,:) = LN - mu*speye(N);
Lap(1+N:5*N,:) = LM - mu*[speye(N); speye(N); speye(N); speye(N)];

%Pre-defining the size of variables, for speed
F = zeros(5*N,1);
Q = zeros(5*N,1);
C = zeros(5*N,1);

% Quadratic nonlinearity
Q(1:N) = u0.^2 + 2*u1.^2 + 2*u2.^2 + 2*u3.^2 + 2*u4.^2;
Q(N+1:2*N) = 2.*(u0.*u1 + u2.*u1 + u2.*u3 + u3.*u4);
Q(2*N+1:3*N) = u1.^2 + 2.*(u2.*u0 + u1.*u3 + u2.*u4);
Q(3*N+1:4*N) = 2.*(u1.*u4 + u1.*u2 + u0.*u3);
Q(4*N+1:5*N) = u2.^2 + 2.*(u1.*u3 + u0.*u4);

% Cubic nonlinearity
C(1:N) = u0.^3 + 6*(u0.*(u1.^2 + u2.^2 + u3.^2 + u4.^2) + u2.*u1.^2 + u4.*u2.^2 + 2.*(u2 + u4).*u1.*u3);
C(N+1:2*N) = (3*(u1.^2 + u0.^2) + 6*(u2.^2 + u2.*u0 + u3.^2 + u4.^2 + u2.*u4)).*u1 + (3.*(u1.^2 + u2.^2) + 6.*(u0.*u2 + u0.*u4 + u2.*u4)).*u3;
C(2*N+1:3*N) = 3*(u2.^3 + u2.*u0.^2 + u0.*u1.^2 + u4.*u1.^2 + u4.*u3.^2) + 6*(u2.*u1.^2 + u2.*u3.^2 + u2.*u4.^2 + u0.*u2.*u4 + u0.*u1.*u3 + u1.*u2.*u3 + u1.*u3.*u4);
C(3*N+1:4*N) = u1.^3 + 3.*(u3.*u0.^2 + u1.*u2.^2 + u3.^3) + 6.*(u0.*u1.*u2 + u0.*u1.*u4 + u3.*u1.^2 + u1.*u2.*u4 + u3.*u2.^2 + u2.*u3.*u4 + u3.*u4.^2);
C(4*N+1:5*N) = 3.*(u4.*u0.^2 + u0.*u2.^2 + u2.*u1.^2 + u2.*u3.^2 + u4.^3) + 6.*(u0.*u1.*u3 + u4.*u1.^2 + u1.*u2.*u3 + u4.*u2.^2 + u4.*u3.^2);

% F will be the final output function
F(1:N) = Lap(1:N,1:N)*u0 + p(2)*Q(1:N) + p(3)*C(1:N);
F(N+1:2*N) = Lap(1+ N:2*N,1:N)*u1 + p(2)*Q(1+N:2*N)+ p(3)*C(1+N:2*N);
F(2*N+1:3*N) = Lap(1+ 2*N:3*N,1:N)*u2 + p(2)*Q(1+2*N:3*N) + p(3)*C(1+2*N:3*N);
F(3*N+1:4*N) = Lap(1+ 3*N:4*N,1:N)*u3 + p(2)*Q(3*N+1:4*N) + p(3)*C(3*N+1:4*N);
F(4*N+1:5*N) = Lap(1+ 4*N:5*N,1:N)*u4 + p(2)*Q(4*N+1:5*N) + p(3)*C(4*N+1:5*N);

if nargout > 1 %Calculating the Jacobian
        
J = zeros(5*N,5*N);
DLap = zeros(5*N,5*N);
DQ = zeros(5*N,5*N);
DC = zeros(5*N,5*N);

DLap(1:N,1:N) = Lap(1:N,1:N);
DLap(N+1:2*N,N+1:2*N) = Lap(N+1:2*N,1:N);
DLap(2*N+1:3*N,2*N+1:3*N) = Lap(2*N+1:3*N,1:N);
DLap(3*N+1:4*N,3*N+1:4*N) = Lap(3*N+1:4*N,1:N);
DLap(4*N+1:5*N,4*N+1:5*N) = Lap(4*N+1:5*N,1:N);

DQ(1:N,1:N) = 2*sparse(1:N,1:N,u0);
DQ(1:N,N+1:2*N) = 4.*sparse(1:N,1:N,u1);
DQ(1:N,2*N+1:3*N) = 4.*sparse(1:N,1:N,u2);
DQ(1:N,3*N+1:4*N) = 4.*sparse(1:N,1:N,u3);
DQ(1:N,4*N+1:5*N) = 4.*sparse(1:N,1:N,u4);

DQ(N+1:2*N,1:N) = 2*sparse(1:N,1:N,u1);
DQ(N+1:2*N,N+1:2*N) = 2.*sparse(1:N,1:N,u0 + u2);
DQ(N+1:2*N,2*N+1:3*N) = 2.*sparse(1:N,1:N,u1 + u3);
DQ(N+1:2*N,3*N+1:4*N) = 2.*sparse(1:N,1:N, u2 + u4);
DQ(N+1:2*N,4*N+1:5*N) = 2.*sparse(1:N,1:N, u3);

DQ(2*N+1:3*N,1:N) = 2*sparse(1:N,1:N,u2);
DQ(2*N+1:3*N,N+1:2*N) = 2.*sparse(1:N,1:N,u1 + u3);
DQ(2*N+1:3*N,2*N+1:3*N) = 2.*sparse(1:N,1:N,u0 + u4);
DQ(2*N+1:3*N,3*N+1:4*N) = 2.*sparse(1:N,1:N,u1);
DQ(2*N+1:3*N,4*N+1:5*N) = 2.*sparse(1:N,1:N,u2);

DQ(3*N+1:4*N,1:N) = 2*sparse(1:N,1:N,u3);
DQ(3*N+1:4*N,N+1:2*N) = 2.*sparse(1:N,1:N,u2 + u4);
DQ(3*N+1:4*N,2*N+1:3*N) = 2.*sparse(1:N,1:N,u1);
DQ(3*N+1:4*N,3*N+1:4*N) = 2.*sparse(1:N,1:N,u0);
DQ(3*N+1:4*N,4*N+1:5*N) = 2.*sparse(1:N,1:N,u1);

DQ(4*N+1:5*N,1:N) = 2*sparse(1:N,1:N,u4);
DQ(4*N+1:5*N,N+1:2*N) = 2.*sparse(1:N,1:N,u3);
DQ(4*N+1:5*N,2*N+1:3*N) = 2.*sparse(1:N,1:N,u2);
DQ(4*N+1:5*N,3*N+1:4*N) = 2.*sparse(1:N,1:N,u1);
DQ(4*N+1:5*N,4*N+1:5*N) = 2.*sparse(1:N,1:N,u0);

DC(1:N,1:N) = 3.*sparse(1:N,1:N, u0.^2 + 2.*u1.^2 + 2.*u2.^2 + 2.*u3.^2 + 2.*u4.^2);
DC(1:N,N+1:2*N) = 12.*sparse(1:N,1:N, u0.*u1 + u1.*u2 + u2.*u3 + u3.*u4);
DC(1:N,2*N+1:3*N) = 6.*sparse(1:N,1:N, 2.*u0.*u2 + u1.^2 + 2.*u1.*u3 + 2.*u2.*u4);
DC(1:N,3*N+1:4*N) = 12.*sparse(1:N,1:N, u0.*u3 + u1.*u2 + u1.*u4);
DC(1:N,4*N+1:5*N) = 6.*sparse(1:N,1:N, u2.^2 + 2.*u0.*u4 + 2.*u1.*u3);

DC(N+1:2*N,1:N) = 6.*sparse(1:N,1:N, u0.*u1 + u1.*u2 + u2.*u3 + u3.*u4);
DC(N+1:2*N,N+1:2*N) = 3.*sparse(1:N,1:N, u0.^2 + 2.*u0.*u2 + 3.*u1.^2 + 2.*u2.^2 + 2.*u1.*u3 + 2.*u2.*u4 + 2.*u3.^2 + 2.*u4.^2);
DC(N+1:2*N,2*N+1:3*N) = 6.*sparse(1:N,1:N, u0.*u1 + 2.*u1.*u2 + u0.*u3 + u1.*u4 + u2.*u3 + u3.*u4);
DC(N+1:2*N,3*N+1:4*N) = 3.*sparse(1:N,1:N, 2.*u0.*u2 + 2.*u0.*u4 + u1.^2 + 4.*u1.*u3 + u2.^2 + 2.*u2.*u4);
DC(N+1:2*N,4*N+1:5*N) = 6.*sparse(1:N,1:N, u0.*u3 + u1.*u2 + 2.*u1.*u4 + u2.*u3);

DC(2*N+1:3*N,1:N) = 3.*sparse(1:N,1:N, 2.*u0.*u2 + u1.^2 + 2.*u1.*u3 + 2.*u2.*u4);
DC(2*N+1:3*N,N+1:2*N) = 6.*sparse(1:N,1:N, u0.*u1 + 2.*u1.*u2 + u0.*u3 + u1.*u4 + u2.*u3 + u3.*u4);
DC(2*N+1:3*N,2*N+1:3*N) = 3.*sparse(1:N,1:N, u0.^2 + 2.*u1.^2 + 3.*u2.^2 + 2.*u0.*u4 + 2.*u1.*u3 + 2.*u3.^2 + 2.*u4.^2);
DC(2*N+1:3*N,3*N+1:4*N) = 6.*sparse(1:N,1:N, u0.*u1 + u1.*u2 + u1.*u4 + 2.*u3.*u2 + u3.*u4);
DC(2*N+1:3*N,4*N+1:5*N) = 3.*sparse(1:N,1:N, 2.*u0.*u2 + u1.^2 + 2.*u1.*u3 + 4.*u4.*u2 + u3.^2);

DC(3*N+1:4*N,1:N) = 6.*sparse(1:N,1:N, u3.*u0 + u1.*u2 + u1.*u4);
DC(3*N+1:4*N,N+1:2*N) = 3.*sparse(1:N,1:N, u1.^2 + u2.^2 + 2.*(u0.*u2 + u0.*u4 + u2.*u4) + 4.*u3.*u1);
DC(3*N+1:4*N,2*N+1:3*N) = 6.*sparse(1:N,1:N,u0.*u1 + u1.*u2 + u1.*u4 + u3.*u4 + 2.*u2.*u3);
DC(3*N+1:4*N,3*N+1:4*N) = 3.*sparse(1:N,1:N, u0.^2 + 2.*(u1.^2 + u2.^2 + u2.*u4 + u4.^2) + 3.*u3.^2);
DC(3*N+1:4*N,4*N+1:5*N) = 6.*sparse(1:N,1:N, u0.*u1 + u1.*u2 + u2.*u3 + 2.*u3.*u4);

DC(4*N+1:5*N,1:N) = 3.*sparse(1:N,1:N, u2.^2 + 2.*u0.*u4 + 2.*u1.*u3);
DC(4*N+1:5*N,N+1:2*N) = 6.*sparse(1:N,1:N, u0.*u3 + u1.*u2 + u2.*u3 + 2.*u1.*u4);
DC(4*N+1:5*N,2*N+1:3*N) = 3.*sparse(1:N,1:N, u1.^2 + u3.^2 + 2.*(u0.*u2 + u1.*u3) + 4.*u2.*u4);
DC(4*N+1:5*N,3*N+1:4*N) = 6.*sparse(1:N,1:N, u0.*u1 + u1.*u2 + u2.*u3 + 2.*u3.*u4);
DC(4*N+1:5*N,4*N+1:5*N) = 3.*sparse(1:N,1:N, u0.^2 + 2.*(u1.^2 + u2.^2 + u3.^2) + 3.*u4.^2);


J = DLap + p(2)*DQ + p(3)*DC;
end

end