% Bivariate polynomial interpolation in chebyshev basis. Takes in a
% function handle and the degree of the desired interpolation polynomial
function result = bivar_cheby_interpolate(fun,n)

% 2D Chebyshev points
[x,y] = ndgrid(cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)');

% Have to perform some shifts since the Chebyshev points do not start at
% cos(0) but cos(pi/(2*n)). Compute the residuals needed for these shifts.
residuals = (0:n-1)'*pi/(2*n);
residuals = [residuals; 0; -flip(residuals(2:end))];
residuals = exp(-1i*residuals)*exp(-1i*residuals');

% Evaluate the function at the given points (x,y)
f_x = fun(x(:), y(:));
% Express the values in a grid instead of a list
f_x = reshape(f_x,n,n);

% Find the Chebyshev coefficients by using multidimensional fft
cheby_coeffs = fftn([f_x flip(f_x,2); flip(f_x,1) flip(flip(f_x,1),2)])/(4*n^2);
% Perform shifts
cheby_coeffs = cheby_coeffs.*residuals;
% Get the cosine coefficients
cheby_coeffs(:,2:n) = cheby_coeffs(:,2:n) + flip(cheby_coeffs(:,n+2:end),2);
cheby_coeffs(2:n,1:n) = cheby_coeffs(2:n,1:n) + flip(cheby_coeffs(n+2:end,1:n),1);

% Matrix representation of a bivariate polynomial in Chebyshev basis.
cheby_coeffs = cheby_coeffs(1:n,1:n);
% Degrees in decreasing order as in the input polynomial
result = flip(flip(cheby_coeffs, 1), 2);

% Code below is for debugging

% Chebyshev polynomials in monomial basis
A = cheby_basis(n);

% Sanity check: Check that we get the the interpolation makes sense in monomial basis
monomial_coeffs = A*cheby_coeffs*A'

end