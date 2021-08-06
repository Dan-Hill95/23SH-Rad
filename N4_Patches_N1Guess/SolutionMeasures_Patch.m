% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function F = SolutionMeasures_Patch(step,uu,pp,mesh_params)

  %% Rename parameters
     r = mesh_params.r;
     N = mesh_params.N;
    Lr = max(r);
    hr = abs(r(2)-r(1));
     u = uu(1:5*N);
  %% Compute quadrature weights and l2norm
     w = ones(N,1); 
l2Norm0 = sum( hr * w .* u(1+0*N:1*N).^2)/(2*Lr);
l2Norm1 = sum( hr * w .* u(1+1*N:2*N).^2)/(2*Lr);
l2Norm2 = sum( hr * w .* u(1+2*N:3*N).^2)/(2*Lr);
l2Norm3 = sum( hr * w .* u(1+3*N:4*N).^2)/(2*Lr);
l2Norm4 = sum( hr * w .* u(1+4*N:5*N).^2)/(2*Lr);

  %% Allocate
     F = zeros(1,5);

  %% Assign branch variables
     F = [l2Norm0 l2Norm1 l2Norm2 l2Norm3 l2Norm4];
  
end
