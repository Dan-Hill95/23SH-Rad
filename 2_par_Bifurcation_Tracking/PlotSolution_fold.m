% Adapted by Dan J Hill 2021
% Original Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function plotHandle = PlotSolution_fold(u,parentHandle,mesh_params)
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
   uu=u(1:N);
   uv=u(1+N:2*N);
  subplot(1,2,1)        %Plot of solution
  plot(r,uu(1:N),'.-');
   pbaspect([1 1 1]);
   axis([0 50 -1.2*max(-uu) 1.2 *max(uu)]);
  subplot(1,2,2)        % Plot of destabilising eigenmode
   plot(r,uv(1:N),'r--');
   pbaspect([1 1 1]);
   axis([0 50 -1.2*max(-uv) 1.2 *max(uv)]);
   
   %% Save
   % print -depsc state.eps
   print -dtiff state.tiff

end

