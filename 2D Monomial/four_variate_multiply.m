% Computes the multiplication of two 4-variate polynomials A and B in
% tensor representation.
function product = four_variate_multiply(A,B)
    [n1, n2, n3, n4] = size(A);
    [m1, m2, m3, m4] = size(B);
    product = zeros(n1+m1-1,n2+m2-1,n3+m3-1,n4+m4-1);
    for i = 1:n4
        for j = 1:m4
            product(:,:,:,i+j-1) = product(:,:,:,i+j-1)+trivariate_multiply(A(:,:,:,i),B(:,:,:,j));
        end
    end
    
end

