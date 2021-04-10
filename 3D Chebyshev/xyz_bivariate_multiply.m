% Computes the multiplication of two bivariate polynomials A and B in
% matrix representation. A has variables x and y, and B has variables z and
% y. Returns a tensor that has variables x,y,z.
function product = xyz_bivariate_multiply(A,B)
    [n1, n2] = size(A);
    [m1, m2] = size(B);
    product = zeros(n1,n2+m2-1,m1);
    for i = 1:n1
        product(i,:,:) = bivariate_multiply(A(i,:),B)';
    end
    
end

