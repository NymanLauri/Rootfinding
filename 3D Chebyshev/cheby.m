% Univariate polynomial interpolation in chebyshev basis

% Vector representation of a polynomial whose Chebyshev coefficients we
% want to find.
poly=[1 + 2*i; 7; 5];

n = length(poly);
x = cos((2*(1:n)-1)/(2*n)*pi)';

% This choice of points also works with small tweaks in the code
% x = cos((0:n-1)/n*pi)'

residuals = (1:n-1)'*pi/(2*n);

f_x = poly_eval(poly, x);

cheby_coeffs = fft([f_x; flip(f_x)])/(2*n);
cheby_coeffs(2:n) = cheby_coeffs(2:n).*exp(-residuals*i) + flip(cheby_coeffs(n+2:end)).*exp(residuals*i);
cheby_coeffs = cheby_coeffs(1:n);

A = cheby_basis(n)

A*cheby_coeffs
