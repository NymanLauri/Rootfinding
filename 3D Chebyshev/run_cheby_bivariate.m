% Bivariate polynomial interpolation in chebyshev basis

% Matrix representation of a bivariate polynomial whose Chebyshev coefficients we
% want to find. Degrees in decreasing order. Assumes that the matrix is
% square.
poly=[0 0 0 0; 0 4 0 -2; 0 0 0 0; 0 -2 0 1];
n = length(poly);

% Evaluate the polynomial at the given points (x,y)
fun = @(x,y) poly_eval_bivariate(poly, [x y]);

% fun = @(x,y) exp(x)
% n = 10

bivar_cheby_interpolate(fun,n)