% 5-variate polynomial interpolation in chebyshev basis. Takes in the
% values of a function at Chebyshev points and returns a tensor
% representing the interpolation polynomial in Chebyshev basis.
function result = cheby_5D_interpolate(f)
[n1,n2,n3,n4,n5] = size(f);

% Using symmetricity of cosine
f = [f flip(f,2); flip(f,1) flip(flip(f,1),2)];
f(:,:,n3+1:2*n3,:,:) = flip(f,3);
f(:,:,:,n4+1:2*n4,:) = flip(f,4);
f(:,:,:,:,n5+1:2*n5) = flip(f,5);

% Have to perform some shifts since the Chebyshev points do not start at
% cos(0) but cos(pi/(2*n)). Compute the residuals needed for these shifts.
residuals1 = (0:n1-1)'*pi/(2*n1);
residuals1 = [residuals1; 0; -flip(residuals1(2:end))];
residuals2 = (0:n2-1)'*pi/(2*n2);
residuals2 = [residuals2; 0; -flip(residuals2(2:end))];
residuals3 = (0:n3-1)'*pi/(2*n3);
residuals3 = [residuals3; 0; -flip(residuals3(2:end))];
residuals4 = (0:n4-1)'*pi/(2*n4);
residuals4 = [residuals4; 0; -flip(residuals4(2:end))];
residuals5 = (0:n5-1)'*pi/(2*n5);
residuals5 = [residuals5; 0; -flip(residuals5(2:end))];

temp = exp(-1i*residuals1)*exp(-1i*residuals2');
temp = kron(exp(-1i*residuals3'),temp);
temp = kron(exp(-1i*residuals4'),temp);
temp = kron(exp(-1i*residuals5'),temp);
residuals = reshape(temp,2*n1,2*n2,2*n3,2*n4,2*n5);

% Find the Chebyshev coefficients by using multidimensional fft
cheby_coeffs = fftn(f)/(32*n1*n2*n3*n4*n5);
% Perform shifts
cheby_coeffs = cheby_coeffs.*residuals;
% Get the cosine coefficients
cheby_coeffs(:,2:n2,:,:,:) = cheby_coeffs(:,2:n2,:,:,:) + flip(cheby_coeffs(:,n2+2:end,:,:,:),2);
cheby_coeffs(2:n1,1:n2,:,:,:) = cheby_coeffs(2:n1,1:n2,:,:,:) + flip(cheby_coeffs(n1+2:end,1:n2,:,:,:),1);
cheby_coeffs(1:n1,1:n2,2:n3,:,:) = cheby_coeffs(1:n1,1:n2,2:n3,:,:) + flip(cheby_coeffs(1:n1,1:n2,n3+2:end,:,:),3);
cheby_coeffs(1:n1,1:n2,1:n3,2:n4,:) = cheby_coeffs(1:n1,1:n2,1:n3,2:n4,:) + flip(cheby_coeffs(1:n1,1:n2,1:n3,n4+2:end,:),4);
cheby_coeffs(1:n1,1:n2,1:n3,1:n4,2:n5) = cheby_coeffs(1:n1,1:n2,1:n3,1:n4,2:n5) + flip(cheby_coeffs(1:n1,1:n2,1:n3,1:n4,n5+2:end),5);

% Tensor representation of a trivariate polynomial in Chebyshev basis.
% Degrees in increasing order.
result = cheby_coeffs(1:n1,1:n2,1:n3,1:n4,1:n5);

end