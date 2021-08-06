% Adapted by Dan J Hill 2021
% from codes by David J.B. Lloyd - University of Surrey

%% fold continuation
function [F,J] = bif_cont_SH(uu,p,mesh_params,problem_handle)
n=mesh_params.N;                                                           % Length of the mesh
LN = mesh_params.LN;                                                       % Linear radial operator
LM = mesh_params.LM;                                                       % Linear D_{m} operator
ContVar=mesh_params.ContVar;                                               % Index of parameter to be varied 
u  = uu(1:n);                                                              % Radial solution
v = uu(n+1:2*n);                                                           % D_{m} perturbation
pp = p;                                                                    % Calling the parameters of the system
pp(ContVar) = uu(2*n+1);                                                   % Updating the varied parameter  
                                                                           %
[Fu,Ju] = problem_handle(uu(1:n),pp);                                      % Compute the Jacobian of the function
                                                                           %
Jv = Ju - LN + LM;                                                         % Replace radial perturbations with D_{m} perturbations
                                                                           %  
F = [Fu; Jv*v; norm(v,2)^2-1];                                             % Append the linear stability problem for normalised eigenmodes
                                                                           %   
if nargout > 1                                                             % Calculate the extended Jacobian (only when called for)
                                                                           %
J = sparse(2*n+1,2*n+1);                                                                           %
Juu = spdiags((2*p(2) + 6*p(3)*u).*v,0,n,n);                               % Second-order variations of F with respect to u   
                                                                           %
J= [ Ju sparse(n,n) zeros(n,1); Juu Jv zeros(n,1); zeros(1,n) 2*v' 0];     %    
                                                                           %
epsiF = 1e-8;                                                              % Perturbing by the parameter p(ContVar)
p2 = pp; p2(ContVar) = p2(ContVar)+epsiF;                                  %
dF = problem_handle(uu(1:n),p2);                                           %
J(1+n:2*n,2*n+1) = 2*u.*v;                                                 %
J(1:n,2*n+1) = (dF - Fu)/epsiF;                                            %
end         