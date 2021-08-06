% Adapted by Dan J Hill 2021
% Original: Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function F = SolutionMeasures_SH(uu,mesh_params)

  %% Rename parameters
     r = mesh_params.r;
     N = mesh_params.N;
    Lr = max(r);
    hr = abs(r(2)-r(1));
     u = uu(1:N);
  %% Compute quadrature weights and l2norm
     w = ones(size(u)); 
l2Norm = sum( hr * w .* u.^2)/(2*Lr);

  %% Allocate
     F = zeros(1,1);

  %% Assign branch variables
     F = [l2Norm];
  
end
