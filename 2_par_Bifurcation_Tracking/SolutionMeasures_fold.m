% Adapted by Dan J Hill 2021
% Original: Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function F = SolutionMeasures_fold(u,mesh_params)

  %% Rename parameters
   x = mesh_params.r;
   N = mesh_params.N;
  Lx = max(x);
  hx = abs(x(2)-x(1));
  uu=u(1:N);
  uv=u(1+N:2*N);

  %% Compute quadrature weights and l2norm
  w = ones(size(uu)); 
  l2Normu = sum( hx * w .* uu.^2)/(2*Lx);
  l2Normv = sum( hx * w .* uv.^2)/(2*Lx);

  %% Allocate
  F = zeros(1,3);

  %% Assign branch variables
  F = [l2Normu l2Normv u(2*N+1)];
  
end
