% Adapted by Dan J Hill 2021
% Original Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function [V,LAMBDA] = ComputeEigenvalues_SH(uu,p,mesh_params)
    LN = mesh_params.LN;    % Calling radial Linear Operator
    LM = mesh_params.LM;    % Calling D_{m} Linear Operator
  %% Compute linear operators
  [~,J] = Equation_SH(uu,p,mesh_params);    % Computing the Jacobian of the Swift-Hohenberg equation at solution uu
  J = J - LN + LM;                          % Replacing radial perturbations with D_{m} perturbations
  
  %% Call direct eigenvalue solver
  [V,LAMBDA] = eigs(full(J),10,1);        % Only compute the 10 nearest eigenvalues to the right half-plane (i.e. with positive real part)
end
