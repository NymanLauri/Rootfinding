% Evaluate the polynomial at the given points. Assumes that coeffs is a
% square matrix.
function y = poly_eval_bivariate(coeffs,points)
n = size(coeffs,1);
n_points = size(points,1);

a = vander(points(:,1));
a = a(:,end-n+1:end);

b = vander(points(:,2));
b = b(:,end-n+1:end);

y = repmat(a,1,n).*reshape(repmat(b,n,1),n_points,n^2)*coeffs(:);
end
