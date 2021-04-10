% Computes the multiplication of two bivariate polynomials A and B in
% matrix representation.
function product = bivariate_multiply(A,B)
    [n1, n2] = size(A);
    [m1, m2] = size(B);
    product = zeros(n1+m1-1,n2+m2-1);
    for i = 1:n2
        for j = 1:m2
            product(:,i+j-1) = product(:,i+j-1)+univariate_multiply(A(:,i),B(:,j));
        end
    end
    
end

