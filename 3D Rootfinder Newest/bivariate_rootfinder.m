% Takes in two bivariate functions whose zeros we want to find.
% Solves the second component of the roots.
% n is the maximal degree of the polynomials.
function roots = bivariate_rootfinder(f1,f2,n)

% The Cayley function
f_cayley = @(s,t,y) (f1(s,y).*f2(t,y) - f1(t,y).*f2(s,y)) ./ (s-t);

% Amounts of interpolation points of s,t,y
n_s = n;
% In order to avoid dividing by 0, add +1
n_t = n+1;
n_y = 2*n+1;

% Interpolation points 
s_vals = cos((2*(1:n_s)-1)/(2*n_s)*pi)';
t_vals = cos((2*(1:n_t)-1)/(2*n_t)*pi)';
y_vals = cos((2*(1:n_y)-1)/(2*n_y)*pi)';

% 5D Chebyshev points
[a,b,c] = ndgrid(s_vals,t_vals,y_vals);

% Evaluate the function at the given points (a,b,c,d,e)
f = f_cayley(a(:), b(:), c(:));
% Express the values in a grid instead of a list
f = reshape(f,n_s,n_t,n_y);

R = cheby_3D_interpolate(f);
% The degree of t is actually n so we can leave column n+1 out
R = R(1:n,1:n,:);

% If all the coefficients of the highest degree are zero, leave it out from the linearization 
y_length = n_y;
for i = flip(1:n_y)
    if (norm(R(:,:,i),'fro') < 1e-12); y_length = i-1; else; break; end
end

if (y_length == 0)
    roots = [];
elseif (y_length <= 2)
    [~,D] = eig(R(:,:,1),-R(:,:,2));
    roots = diag(D);
else
    % Compute the linearization matrices C1 and C2
    C1 = -2*eye(n*(y_length-1));
    C1(end-n+1:end,end-n+1:end) = 2*R(:,:,y_length);
    C1(1:n,1:n) = -eye(n);

    C2 = zeros(size(C1));
    C2(1:end-n,n+1:end) = eye(n*(y_length-2));
    C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)+eye(n*(y_length-2));
    C2(end-n+1:end,end-2*n+1:end-n) = -R(:,:,y_length);
    
    % Compute the last rows of the coefficient matrix C2
    D=[];
    for i = 1:y_length-1
        D = [D R(:,:,i)];
    end
    C2(end-n+1:end,:) = C2(end-n+1:end,:)+D;

    % Solve the eigenproblem
    [~,D] = eig(C2,-C1);
    
%     % Solve the eigenproblem
%     [~,D] = eig(vpa(C2),-vpa(C1));

    roots = diag(D);
    % Ignore the roots outside the interval [-1,1] (or the complex unit disk at
    % this point)
    roots = roots(abs(roots) < 1);
end
end
