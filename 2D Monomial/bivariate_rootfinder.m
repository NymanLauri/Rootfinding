% Computes the first components of the roots of the bivariate polynomials 
% p and q that are in matrix form. The algorithm finds the roots of the variable
% whose degree is determined by the column number in the matrix representation.
function [V,D] = bivariate_rootfinder(p,q)
    % Assumes that p and q are the same size and square.
    assert(all(size(p)==size(q)))
    assert(size(p,1)==size(p,2))

    n = size(p,1);

    % Compute the Bezout matrix B
    B_initial = xyz_bivariate_multiply(p,q) - xyz_bivariate_multiply(q,p);
    % Divide by x-z
    B_final = xz_divide(B_initial);

    % Note: B possibly is too long in the y-direction (coefficients of highest
    % degrees are all zeros). This does not seem to affect the computed
    % eigenvalues.

    B = permute(B_final,[1,3,2]);

    % Compute the linearization matrices C1 and C2
    C1 = eye((n-1)*(size(B,3)-1));
    C1(end-n+2:end,end-n+2:end) = B(:,:,end);

    C2 = zeros(size(C1));
    C2(1:end-n+1,n:end) = -eye((n-1)*(size(B,3)-2));

    % Compute the last rows of the coefficient matrix C2
    D=[];
    for i = 1:size(B,3)-1
        D = [D B(:,:,i)];
    end
    C2(end-n+2:end,:) = D;

    % Solve the eigenproblem
    [V,D] = eig(C2,-C1);
end