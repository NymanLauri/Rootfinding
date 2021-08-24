% Takes in three trivariate functions whose zeros we want to find.
% Solves the third component of the roots in vpa (requires tweaking the 
% code - never got this to work properly, see exact_resultant.m).
% n is the maximal degree of the polynomials.
function roots = trivariate_rootfinder(f1,f2,f3,n)
% Perturbation
f1 = @(x,y,z) f1(x,y,z) + eps*x.^n + eps*y.^n + eps*z.^n;
f2 = @(x,y,z) f2(x,y,z) + eps*x.^n + eps*y.^n + eps*z.^n;
f3 = @(x,y,z) f3(x,y,z) + eps*x.^n + eps*y.^n + eps*z.^n;

% addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
% f1 = chebfun3(@(x,y,z) f1(x,y,z), [n+1, n+1, n+1]);
% f2 = chebfun3(@(x,y,z) f2(x,y,z), [n+1, n+1, n+1]);
% f3 = chebfun3(@(x,y,z) f3(x,y,z), [n+1, n+1, n+1]);

disp('Evaluating the Cayley function:')
tic
% The Cayley function
f_cayley = @(s1,s2,t1,t2,z) (f1(s1,s2,z).*f2(t1,s2,z).*f3(t1,t2,z) + f2(s1,s2,z).*f3(t1,s2,z).*f1(t1,t2,z) + f3(s1,s2,z).*f1(t1,s2,z).*f2(t1,t2,z) ...
    - f3(s1,s2,z).*f2(t1,s2,z).*f1(t1,t2,z) - f2(s1,s2,z).*f1(t1,s2,z).*f3(t1,t2,z) - f1(s1,s2,z).*f3(t1,s2,z).*f2(t1,t2,z)) ./ ((s1-t1).*(s2-t2));

% Amounts of interpolation points based on degrees of s1,s2,t1,t2,z
n_s1 = n;
n_s2 = 2*n;
n_t1 = 2*n;
n_t2 = n;
n_z = 3*n+1;

% Interpolation points 
s1_vals = cos((2*(1:n_s1)-1)/(2*n_s1)*pi)';
s2_vals = cos((2*(1:n_s2)-1)/(2*n_s2)*pi)';
t1_vals = s2_vals;
t2_vals = s1_vals;
z_vals = cos((2*(1:n_z)-1)/(2*n_z)*pi)';

% 5D Chebyshev points
[a,b,c,d,e] = ndgrid(s1_vals,s2_vals,t1_vals,t2_vals,z_vals);

% Evaluate the function at the given points (a,b,c,d,e)
f = f_cayley(a(:), b(:), c(:), d(:), e(:));
% Express the values in a grid instead of a list
f = reshape(f,n_s1,n_s2,n_t1,n_t2,n_z);
toc

disp('Interpolating the Cayley function:')
tic
A = cheby_5D_interpolate(f);
toc

disp('Solving the eigenproblem:')
tic
% Matricization of the tensor
R = reshape(A, n_s1*n_s2, n_s1*n_s2, n_z);

%0.5*(s1+t1+u)*(2 s2 z^2 + 2 t2 z^2 - 2 u z^2 + sqrt(2) u z^2 + sqrt(2) s2 t2 u - sqrt(2) s2 u z - sqrt(2) t2 u z)

% If all the coefficients of the highest degree are zero, leave it out from the linearization 
z_length = n_z;
for i = flip(1:n_z)
    if (norm(R(:,:,i),'fro') < 1e-12); z_length = i-1; else; break; end
end

if (z_length == 0)
    roots = [];
elseif (z_length <= 2)
    [~,D] = eig(R(:,:,1),-R(:,:,2));
    roots = diag(D);
else
    n = n_s1*n_s2;

    % Compute the linearization matrices C1 and C2
    C1 = -2*eye(n*(z_length-1));
    C1(end-n+1:end,end-n+1:end) = 2*R(:,:,z_length);
    C1(1:n,1:n) = -eye(n);

    C2 = zeros(size(C1));
    C2(1:end-n,n+1:end) = eye(n*(z_length-2));
    C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)+eye(n*(z_length-2));
    C2(end-n+1:end,end-2*n+1:end-n) = -R(:,:,z_length);
    
    % Compute the last rows of the coefficient matrix C2
    D=[];
    for i = 1:z_length-1
        D = [D R(:,:,i)];
    end
    C2(end-n+1:end,:) = C2(end-n+1:end,:)+D;

    % Solve the eigenproblem
    [~,D] = eig(C2,-C1);
%     [~, D, C]=polyeig(C2,-C1);

%     C2=vpa(C2,16);
%     C1=vpa(C1,16);
% 
%     [~,D] = eig(-C2/C1);
%     D = double(D);
%     
    roots = diag(D);
    % Ignore the roots outside the interval [-1,1] (or the complex unit disk at
    % this point)
    roots = roots(abs(roots) < 1);
    
end
toc


% tic
% cheby_vectors = flip(cheby_basis(z_length),1);
% temp=roots;
% roots=[];
% for j=1:size(temp)
%     z = temp(j);
%     z_vals = [1];
%     for i=2:z_length
%         z_vals = [z_vals z^(i-1)];
%     end
% 
%     M = R(:,:,1);
%     for i=2:z_length
%         M = M + R(:,:,i)*(z_vals*cheby_vectors(:,i));
%     end
%     
%     det(M)
%     if (det(M) < 1e-6); roots = [roots; z]; end
% end
% toc
end
