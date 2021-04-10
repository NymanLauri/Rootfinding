% Trivariate polynomial interpolation in chebyshev basis. Takes in a
% function handle and the degree of the desired interpolation polynomial
function result = trivar_cheby_interpolate(fun,n)

% 3D Chebyshev points
[x,y,z] = ndgrid(cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)');

% Have to perform some shifts since the Chebyshev points do not start at
% cos(0) but cos(pi/(2*n)). Compute the residuals needed for these shifts.
residuals = (0:n-1)'*pi/(2*n);
residuals = [residuals; 0; -flip(residuals(2:end))];
temp = exp(-1i*residuals)*exp(-1i*residuals');
residuals = kron(exp(-1i*residuals),temp);
residuals = reshape(residuals,2*n,2*n,2*n);

% Evaluate the function at the given points (x,y,z)
f_x = fun(x(:), y(:), z(:));
% Express the values in a grid instead of a list
f_x = reshape(f_x,n,n,n);

% Using symmetricity
f_x = [f_x flip(f_x,2); flip(f_x,1) flip(flip(f_x,1),2)];
f_x(:,:,n+1:2*n) = flip(f_x,3);

% Find the Chebyshev coefficients by using multidimensional fft
cheby_coeffs = fftn(f_x)/(8*n^3);
% Perform shifts
cheby_coeffs = cheby_coeffs.*residuals;
% Get the cosine coefficients
cheby_coeffs(:,2:n,:) = cheby_coeffs(:,2:n,:) + flip(cheby_coeffs(:,n+2:end,:),2);
cheby_coeffs(2:n,1:n,:) = cheby_coeffs(2:n,1:n,:) + flip(cheby_coeffs(n+2:end,1:n,:),1);
cheby_coeffs(1:n,1:n,2:n) = cheby_coeffs(1:n,1:n,2:n) + flip(cheby_coeffs(1:n,1:n,n+2:end),3);

% Tensor representation of a trivariate polynomial in Chebyshev basis.
result = cheby_coeffs(1:n,1:n,1:n);

end