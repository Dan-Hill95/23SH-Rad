% Dan J Hill (2021) - Solving the simple D_{2k} patch matching condition

% In order to run this code, you only need to call the function 'MatchSoln'

% Example input: 

% MatchSoln(0.2*[1,1,1,1,1,1],3,20,0.001);

%% Input: MatchSoln(x,k,r_max,mu)

%     x is the initial guess: an (N+1)-dimensional vector of the form x = [x(1), x(2), ..., x(N+1)]
%     k is the lattice index: solutions are invariant under rotations of pi/k
% r_max is the maximum range of the mesh
%    mu is the localisation of the solution

% For 3|k, a good initial guess will be of the form x=y*[1,...,1], for some small y

% For 3~|k, a good initial guess will be of the form x=y*[-0.5, 1, 1, -0.5,1,1,...,-0.5,1,1], for some small y