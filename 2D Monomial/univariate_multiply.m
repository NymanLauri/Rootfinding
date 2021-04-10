% Computes the multiplication of two univariate polynomials v and w in vector
% representation.
% v and w are assumed to be column vectors
function product = univariate_multiply(v,w)
    n = size(v,1);
    m = size(w,1);
    c = zeros(n+m-1,1);
    c(1:n) = v;
    r = v(1)*eye(1,length(w));
    A = toeplitz(c,r);
    product = A*w;

end
