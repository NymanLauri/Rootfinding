% Takes in three trivariate functions whose zeros we want to find.
% Solves the third component of the roots.
% n is the maximal degree of the polynomials.
function roots = trivariate_rootfinder(f1,f2,f3,n)
% addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')

disp('Evaluating the Cayley function:')
tic
% cheb1 = chebfun3(@(x,y,z) f1(x,y,z));
% cheb2 = chebfun3(@(x,y,z) f2(x,y,z));
% cheb3 = chebfun3(@(x,y,z) f3(x,y,z));
% 
% Dx_cheb1 = diffx(cheb1);
% Dx_cheb2 = diffx(cheb2);
% Dx_cheb3 = diffx(cheb3);
% Dx_cheb = [Dx_cheb1; Dx_cheb2; Dx_cheb3];
% 
% Dy_cheb1 = diffy(cheb1);
% Dy_cheb2 = diffy(cheb2);
% Dy_cheb3 = diffy(cheb3);
% Dy_cheb = [Dy_cheb1; Dy_cheb2; Dy_cheb3];
% 
% Dxy_cheb1 = diffy(diffx(cheb1));
% Dxy_cheb2 = diffy(diffx(cheb2));
% Dxy_cheb3 = diffy(diffx(cheb3));
% Dxy_cheb = [Dxy_cheb1; Dxy_cheb2; Dxy_cheb3];

% The Cayley function
f_cayley = @(s1,s2,t1,t2,z) (f1(s1,s2,z).*f2(t1,s2,z).*f3(t1,t2,z) + f2(s1,s2,z).*f3(t1,s2,z).*f1(t1,t2,z) + f3(s1,s2,z).*f1(t1,s2,z).*f2(t1,t2,z) ...
    - f3(s1,s2,z).*f2(t1,s2,z).*f1(t1,t2,z) - f2(s1,s2,z).*f1(t1,s2,z).*f3(t1,t2,z) - f1(s1,s2,z).*f3(t1,s2,z).*f2(t1,t2,z)) ./ ((s1-t1).*(s2-t2));

% D1_adj = @(s1,s2,t1,t2,z) [0 -(Dx_cheb2(s1,s2,z).*f3(t1,t2,z) - Dx_cheb3(s1,s2,z).*f2(t1,t2,z)) (Dx_cheb2(s1,s2,z).*f3(t1,s2,z) - Dx_cheb3(s1,s2,z).*f2(t1,s2,z)); ...
%     0 (Dx_cheb1(s1,s2,z).*f3(t1,t2,z) - Dx_cheb3(s1,s2,z).*f1(t1,t2,z)) -(Dx_cheb1(s1,s2,z).*f3(t1,s2,z) - Dx_cheb3(s1,s2,z).*f1(t1,s2,z)); ... 
%     0 -(Dx_cheb1(s1,s2,z).*f2(t1,t2,z) - Dx_cheb2(s1,s2,z).*f1(t1,t2,z)) (Dx_cheb1(s1,s2,z).*f2(t1,s2,z) - Dx_cheb2(s1,s2,z).*f1(t1,s2,z))];
% 
% D_f_cayley = @(s1,s2,t1,t2,z,D_M) trace(adjoint([f1(s1,s2,z) f2(s1,s2,z) f3(s1,s2,z); f1(t1,s2,z) f2(t1,s2,z) f3(t1,s2,z); f1(t1,t2,z) f2(t1,t2,z) f3(t1,t2,z)]) ... 
%    * D_M);
% 
% D12_f_cayley = @(s1,s2,t1,t2,z,D2_M,D12_M) trace(adjoint([f1(s1,s2,z) f2(s1,s2,z) f3(s1,s2,z); f1(t1,s2,z) f2(t1,s2,z) f3(t1,s2,z); f1(t1,t2,z) f2(t1,t2,z) f3(t1,t2,z)]) ... 
%     * D12_M + D1_adj(s1,s2,t1,t2,z) * D2_M);

% Amounts of interpolation points based on degrees for s1,s2,t1,t2,z
n_s1 = n;
n_s2 = 2*n;
n_t1 = 2*n;
n_t2 = n;
n_z = 3*n+1;

% The amounts of points need to be the same for computing the derivatives
% n_s1 = 2*n;
% n_s2 = 2*n;
% n_t1 = 2*n;
% n_t2 = 2*n;
% n_z = 3*n+1;

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

% disp('Evaluating the derivatives of the Cayley function matrix:')
% tic
% D1_vals = zeros(n_s1,n_s2,n_t1,n_t2,n_z,3,3);
% D2_vals = zeros(n_s1,n_s2,n_t1,n_t2,n_z,3,3);
% D12_vals = zeros(n_s1,n_s2,n_t1,n_t2,n_z,3,3);
% % Use chebfun to find the values of the derivative?
% for j=1:3
%     f_1 = Dx_cheb(j);
%     f_2 = Dy_cheb(j);
%     f_3 = Dxy_cheb(j);
% 
%     temp1 = f_1(a(:), b(:), e(:));
%     temp2 = f_2(a(:), b(:), e(:));
%     temp3 = f_3(a(:), b(:), e(:));
% 
%     % Express the values in a grid instead of a list
%     temp1 = reshape(temp1,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp2 = reshape(temp2,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp3 = reshape(temp3,n_s1,n_s2,n_t1,n_t2,n_z);
%     
%     D1_vals(:,:,:,:,:,1,j) = temp1;
%     D2_vals(:,:,:,:,:,1,j) = temp2;
%     D12_vals(:,:,:,:,:,1,j) = temp3;
%     
%     
%     temp1 = f_1(c(:), b(:), e(:));
%     temp2 = f_2(c(:), b(:), e(:));
%     temp3 = f_3(c(:), b(:), e(:));
% 
%     % Express the values in a grid instead of a list
%     temp1 = reshape(temp1,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp2 = reshape(temp2,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp3 = reshape(temp3,n_s1,n_s2,n_t1,n_t2,n_z);
%     
%     D1_vals(:,:,:,:,:,2,j) = temp1;
%     D2_vals(:,:,:,:,:,2,j) = temp2;
%     D12_vals(:,:,:,:,:,2,j) = temp3;
%     
%     
%     temp1 = f_1(c(:), d(:), e(:));
%     temp2 = f_2(c(:), d(:), e(:));
%     temp3 = f_3(c(:), d(:), e(:));
% 
%     % Express the values in a grid instead of a list
%     temp1 = reshape(temp1,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp2 = reshape(temp2,n_s1,n_s2,n_t1,n_t2,n_z);
%     temp3 = reshape(temp3,n_s1,n_s2,n_t1,n_t2,n_z);
%     
%     D1_vals(:,:,:,:,:,3,j) = temp1;
%     D2_vals(:,:,:,:,:,3,j) = temp2;
%     D12_vals(:,:,:,:,:,3,j) = temp3;
%     
% end
% toc 
% 
% 
% 
% disp('Evaluating the derivatives of the Cayley function:')
% tic
% for s1=1:n_s1
%    for s2=1:n_s2
%        for z=1:n_z
%            for t1=1:n_t1
% %                keyboard
%                f(s1,s2,t1,s2,z) = D_f_cayley(s1_vals(s1),s2_vals(s2),t1_vals(t1),t2_vals(s2),z_vals(z),reshape(D2_vals(s1,s2,t1,s2,z,:,:),3,3))/(s1_vals(s1)-t1_vals(t1));
%            end
%            for t2=1:n_t2
%                f(s1,s2,s1,t2,z) = D_f_cayley(s1_vals(s1),s2_vals(s2),t1_vals(s1),t2_vals(t2),z_vals(z),reshape(D1_vals(s1,s2,s1,t2,z,:,:),3,3))/(s2_vals(s2)-t2_vals(t2));
%            end
%            f(s1,s2,s1,s2,z) = D12_f_cayley(s1_vals(s1),s2_vals(s2),t1_vals(s1),t2_vals(s2),z_vals(z),reshape(D2_vals(s1,s2,s1,s2,z,:,:),3,3),reshape(D12_vals(s1,s2,s1,s2,z,:,:),3,3));
%        end
%    end
% end
% toc

disp('Interpolating the Cayley function:')
tic
A = cheby_5D_interpolate(f);
toc

% In ideal case
% N = factorial(d-1)*n^(d-1); R = reshape(A, N, N);

disp('Solving the eigenproblem:')
tic
% Matricization of the tensor
R = reshape(A, n_s1*n_s2, n_s1*n_s2, n_z);

% If all the coefficients of the highest degree are zero, leave it out from the linearization 
z_length = n_z;
for i = flip(1:n_z)
    if (norm(R(:,:,i),'fro') < 1e-12); z_length = i-1; else; break; end
end

if (z_length <= 2)
    [~,D] = eig(R(:,:,1),-R(:,:,2));
    roots = diag(D);
else
    n = n_s1*n_s2;

    % Compute the linearization matrices C1 and C2
    % C1 = eye(n*(z_length-1));
    % C1(1:n,1:n) = R(:,:,z_length);
    % 
    % C2 = zeros(size(C1));
    % C2(1:end-n,n+1:end) = -eye(n*(z_length-2));
    % C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)-eye(n*(z_length-2));
    % C2(end-n+1:end, end-2*n+1:end-n) = -2*eye(n,n);
    % 
    % % Compute the first rows of the coefficient matrix C2
    % D=[];
    % for i = flip(1:z_length-1)
    %     D = [D R(:,:,i)];
    % end
    % C2(1:n,:) = C2(1:n,:)+D;
    % C2 = 0.5*C2;

    % Compute the linearization matrices C1 and C2
    C1 = -2*eye(n*(z_length-1));
    C1(end-n+1:end,end-n+1:end) = 2*R(:,:,z_length);
    C1(1:n,1:n) = -eye(n);

    C2 = zeros(size(C1));
    C2(1:end-n,n+1:end) = eye(n*(z_length-2));
    C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)+eye(n*(z_length-2));
    C2(end-n+1:end,end-2*n+1:end-n) = -R(:,:,z_length);
    % C2(end-n+1:end, end-2*n+1:end-n) = -2*eye(n,n);

    % Compute the last rows of the coefficient matrix C2
    D=[];
    for i = 1:z_length-1
        D = [D R(:,:,i)];
    end
    C2(end-n+1:end,:) = C2(end-n+1:end,:)+D;

    % Solve the eigenproblem
    [~,D] = eig(C2,-C1);

    roots = diag(D);
end
toc
end
