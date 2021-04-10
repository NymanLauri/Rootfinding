% Finding the monomial coefficients
poly=[1 0 3 41]';

n = length(poly);
x = [1; i; -1; -i]

f_x = poly_eval(poly,x)

mono_coeffs = fft(f_x)/n

% Matrix representation of a bivariate polynomial whose Chebyshev coefficients we
% want to find. Degrees in decreasing order
poly=[1 1 2i; 7 2i 1; 0 0 3];
n = size(poly,1)

[a,b] = ndgrid(exp(1i*2*pi*(0:n-1)/n),exp(1i*2*pi*(0:n-1)/n));
x = [a(:),b(:)]

f_x = poly_eval_bivariate(poly,x)
f_x = reshape(f_x,n,n)

% Gives the coefficients smallest degree first (inverted order)
mono_coeffs = fftn(f_x)/n^2