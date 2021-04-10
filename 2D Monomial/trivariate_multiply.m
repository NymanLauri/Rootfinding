% Computes the multiplication of two trivariate polynomials A and B in
% tensor representation.
function product = trivariate_multiply(A,B)
    [n1, n2, n3] = size(A);
    [m1, m2, m3] = size(B);
    product = zeros(n1+m1-1,n2+m2-1,n3+m3-1);
    for i = 1:n3
        for j = 1:m3
            product(:,:,i+j-1) = product(:,:,i+j-1)+bivariate_multiply(A(:,:,i),B(:,:,j));
        end
    end
    
end

