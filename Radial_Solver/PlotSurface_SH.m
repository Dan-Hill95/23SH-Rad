% Adapted by Dan J Hill 2021
% Original Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function plotHandle = PlotSurface_SH(u,parentHandle,mesh_params)
r = mesh_params.r;
N = mesh_params.N;

   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);

   %% 1D plot
   
%    plot(r,u(1:N),'.-');
   
   %% 2D plot
   
   t = 0:0.01:2*pi;
   t = t';
  u = u(1:N);
  [R,T]=meshgrid(r,t);
  [U,T0]=meshgrid(u,t);
  subplot(1,2,1) % Contour plot
  z=surf(R.*cos(T), R.*sin(T),U);
  view(0,90);
  z.FaceColor='flat';
  z.FaceAlpha=1;
  z.EdgeColor='none';
   pbaspect([1 1 1]);
  axis([-50*sqrt(2) 50*sqrt(2) -50*sqrt(2) 50*sqrt(2) -1 2]);
  subplot(1,2,2) % Radial profile
   plot(r,u(1:N),'.-');
   pbaspect([1 1 1]);
   axis([0 50 -1.2*max(-u) 1.2 *max(u)]);
   
   %% Save
   % print -depsc state.eps
   print -dtiff state.tiff

end

