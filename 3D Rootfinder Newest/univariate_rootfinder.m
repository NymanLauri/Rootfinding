% Takes in a univariate function whose zeros we want to find.
% n is the maximal degree of the polynomial.
function roots = univariate_rootfinder(f,n)

f_vec = cheby_1D_interpolate(f,n);

z_length = n+1;

% Compute the linearization matrices C1 and C2
C1 = -2*eye(z_length-1);
C1(end,end) = 2*f_vec(z_length);
C1(1,1) = -1;

C2 = zeros(size(C1));
C2(1:end-1,2:end) = eye(z_length-2);
C2(2:end,1:end-1) = C2(2:end,1:end-1)+eye(z_length-2);
C2(end,end-1) = - f_vec(z_length);

% Compute the last row of the coefficient matrix C2
D=[];
for i = 1:z_length-1
    D = [D f_vec(i)];
end
C2(end,:) = C2(end,:)+D;

% Solve the eigenproblem
[~,D] = eig(C2,-C1);

roots = diag(D);

end

