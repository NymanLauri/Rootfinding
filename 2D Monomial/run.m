% Express bivariate polynomials p and q in matrix form s.t. they are the same size and square
% (add zeros if needed to reach the correct shape)

% p=y+2x is given in matrix form as [0 1; 2 0]
% The algorithm finds the y-components of the roots.

% Examples of polynomials whose roots we are interested in:

% n = 4;
% p = rand(n,n);
% q = rand(n,n);

% p = [0 -1 0 0; 4 0 0 0; -6 0 0 0; 3 0 0 0];
% q = [-1 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0];

p = [-1 0 1; 0 0 0; 1 0 0];
q = [0 -1 0; 1 0 0; 0 0 0];

% p = [-1 0 1; 0 0 0; 1 0 0];
% q = [0 -2 1; 0 0 0; 1 0 0];

% p = [-1 0 1 0; 0 0 0 0; 1 0 0 0; 0 0 0 0];
% q = [0 -1 0 0; 1 0 0 0; 0 0 0 0 ; 0 0 0 0];

[V,D] = bivariate_rootfinder(p,q);
roots = diag(D)