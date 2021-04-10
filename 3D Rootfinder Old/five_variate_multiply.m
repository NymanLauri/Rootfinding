% Computes the multiplication of two 5-variate polynomials A and B in
% tensor representation.
function product = five_variate_multiply(A,B)
    [n1, n2, n3, n4, n5] = size(A);
    [m1, m2, m3, m4, m5] = size(B);
    product = zeros(n1+m1-1,n2+m2-1,n3+m3-1,n4+m4-1,n5+m5-1);
    for i = 1:n5
        for j = 1:m5
            product(:,:,:,:,i+j-1) = product(:,:,:,:,i+j-1)+four_variate_multiply(A(:,:,:,:,i),B(:,:,:,:,j));
        end
    end
    
end
