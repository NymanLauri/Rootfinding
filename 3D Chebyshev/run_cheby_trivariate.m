% Trivariate polynomial interpolation in chebyshev basis

% Trivariate function that we want to approximate in Chebyshev basis
fun = @(x,y,z) (2*x.^2 - 1).*(2*y.^2 - 1).*(2*z.^2 - 1) + z - 2*x + 9;
n = 3;

% fun = @(x,y) exp(x)
% n = 10

trivar_cheby_interpolate(fun,n)
