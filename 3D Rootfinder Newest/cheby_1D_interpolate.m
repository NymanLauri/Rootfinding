% Univariate polynomial interpolation in chebyshev basis. Takes in a
% function handle and the degree of the desired interpolation polynomial
function result = cheby_1D_interpolate(fun,n)

n=n+1;
% 2D Chebyshev points
p = cos((2*(1:n)-1)/(2*n)*pi)';

% Have to perform some shifts since the Chebyshev points do not start at
% cos(0) but cos(pi/(2*n)). Compute the residuals needed for these shifts.
residuals = (0:n-1)'*pi/(2*n);
residuals = [residuals; 0; -flip(residuals(2:end))];
residuals = exp(-1i*residuals);

% Evaluate the function at the given points
f_x = fun(p);

% Find the Chebyshev coefficients by using fft
cheby_coeffs = fft([f_x; flip(f_x,1)])/(2*n);
% Perform shifts
cheby_coeffs = cheby_coeffs.*residuals;

% Get the cosine coefficients
cheby_coeffs(2:n) = cheby_coeffs(2:n) + flip(cheby_coeffs(n+2:end),1);

% Vector representation of a univariate polynomial in Chebyshev basis.
cheby_coeffs = cheby_coeffs(1:n);

result = cheby_coeffs;

end