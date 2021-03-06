%Dan J Hill (2021) - Solving the simple D_{2k} patch matching condition
%x is the initial guess: an (N+1)-dimensional vector of the form x = [x(1), x(2), ..., x(N+1)]
%k is the lattice index: solutions are invariant under rotations of pi/k
%r_max is the maximum range of the mesh
%mu is the localisation of the solution

%For 3|k, a good initial guess will be of the form x=y*[1,...,1], for some
%small y

%For 3~|k, a good initial guess will be of the form x=y*[-0.5, 1, 1, -0.5,1,1,...,-0.5,1,1], for some
%small y
function a_out=MatchSoln(x,k,r_max,mu)
%Introduce parameters
N=length(x)-1;
a0=x.*ones(1,N+1);
% r_max=100;
L = 1000;
r = (0:L-1)'*(r_max/(L-1));
t=0:0.01:2*pi;

%fsolve options
options = optimset('Display','iter','TolFun',1e-7);

%Solving the (N+1) D_{2k} matching condition
a_out=fsolve(@(a) match(a,k),a0,options);

%% Plotting the solution: 
%1. Defining the radial amplitudes u[j](r)
for j=1:N+1
u(1+(j-1)*L:j*L) = a_out(j)*(-1)^(j-1)*besselj(2*k*(j-1),r).*exp(-sqrt(mu)*r);
end
%2. Introduce meshes for polar approximation
  [R,T]=meshgrid(r,t);
  [U,T0]=meshgrid(u,t);
%3. Define u(r,t) as polar-Fourier series of radial amplitudes
  UU(:,:,1)=U(:,1:L)/2;
  for j=2:N+1
      UU(:,:,j)=U(:,1+(j-1)*L:j*L);
      UU(:,:,j)=UU(:,:,j).*cos(2*k*(j-1).*T);
  end
%4. Plot contour figure
scrsz = get(0,'ScreenSize');
  z1 = figure('Position',[scrsz(3)/6 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2]);
  k1=surf(R.*cos(T), R.*sin(T),sum(UU,3));
  pbaspect([1 1 1]); 
  view(0,90);
  k1.FaceColor='flat';
  k1.EdgeColor='none';
  set(gca,'Color','none','XColor','none','YColor','none','ZColor','none');
  grid off
  z2 = figure('Position',[scrsz(3)/2 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2]);
  k2=surf(R.*cos(T), R.*sin(T),sum(UU,3));
  pbaspect([1 1 1]); 
  view(45,60);
   k2.FaceColor='flat';
   k2.EdgeColor='none';
   set(gca,'Color','none','XColor','none','YColor','none','ZColor','none');
   grid off
end