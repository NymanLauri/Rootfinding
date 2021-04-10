% Computes the multiplication of two trivariate polynomials A and B in
% tensor representation. A has variables x1, x2 and y, and B has variables z1, z2 and
% y. Returns a tensor that has variables x1,x2,y,z1,z2.
function product = xyz_trivariate_multiply(A,B)
    [n1, n2] = size(A);
    [m1, m2] = size(B);
    product = zeros(n1,n2+m2-1,m1);
    for i = 1:n1
        product(i,:,:) = bivariate_multiply(A(i,:),B)';
    end
    
end

