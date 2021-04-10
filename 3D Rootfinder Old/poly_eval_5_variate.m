% Evaluate the polynomial at the given points. Assumes that coeffs is a
% square tensor.
function y = poly_eval_5_variate(coeffs,points)
n = size(coeffs,1);
n_points = size(points,1);

list_coeffs = flip(coeffs(:),1);

a = vander(points(:,1));
a = a(:,end-n+1:end);

b = vander(points(:,2));
b = b(:,end-n+1:end);

c = vander(points(:,3));
c = c(:,end-n+1:end);

d = vander(points(:,4));
d = d(:,end-n+1:end);

e = vander(points(:,5));
e = e(:,end-n+1:end);

y = repmat(a,1,n^4).*repmat(reshape(repmat(b,n,1),n_points,n^2),1,n^3).*repmat(reshape(repmat(c,n^2,1),n_points,n^3),1,n^2).*repmat(reshape(repmat(d,n^3,1),n_points,n^4),1,n).*reshape(repmat(e,n^4,1),n_points,n^5)*list_coeffs;
end
