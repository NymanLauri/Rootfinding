addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
% addpath('C:\Users\Lauri\Google Drive\Dippa\chebfun-master')

% chebfun('a')

% keyboard
tic
% Three functions in three variables
f1 = @(x,y,z) x+y.*z-1;
f2 = @(x,y,z) z.*(x-y);
f3 = @(x,y,z) (z.^3-0.123456789^3);

n = 5;

% f1 = @(x,y,z) z.*(x-z);
% f2 = @(x,y,z) 2*z.*y.*x;
% f3 = @(x,y,z) z.*(x.^2+y.^2+z.^2-1);
% 
% n=4;

cheb1 = chebfun3(@(x,y,z) f1(x,y,z));
cheb2 = chebfun3(@(x,y,z) f2(x,y,z));
cheb3 = chebfun3(@(x,y,z) f3(x,y,z));

Dx_cheb1 = diffx(cheb1);
Dx_cheb2 = diffx(cheb2);
Dx_cheb3 = diffx(cheb3);
Dx_cheb = [Dx_cheb1; Dx_cheb2; Dx_cheb3];

Dy_cheb1 = diffy(cheb1);
Dy_cheb2 = diffy(cheb2);
Dy_cheb3 = diffy(cheb3);
Dy_cheb = [Dy_cheb1; Dy_cheb2; Dy_cheb3];

Dxy_cheb1 = diffy(diffx(cheb1));
Dxy_cheb2 = diffy(diffx(cheb2));
Dxy_cheb3 = diffy(diffx(cheb3));
Dxy_cheb = [Dxy_cheb1; Dxy_cheb2; Dxy_cheb3];

% p1 = trivar_cheby_interpolate(f1,n);
% p2 = trivar_cheby_interpolate(f2,n);
% p3 = trivar_cheby_interpolate(f3,n);

% A = flip(cheby_basis(n),1);
% 
% temp = reshape(reshape(p1,[],n)*A',n,n,n);
% temp = permute(reshape(reshape(permute(temp,[1 3 2]),[],n)*A',n,n,n),[1 3 2]);
% p1_monomial = permute(reshape(reshape(permute(temp,[3 2 1]),[],n)*A',n,n,n),[3 2 1]);
% 
% temp = reshape(reshape(p2,[],n)*A',n,n,n);
% temp = permute(reshape(reshape(permute(temp,[1 3 2]),[],n)*A',n,n,n),[1 3 2]);
% p2_monomial = permute(reshape(reshape(permute(temp,[3 2 1]),[],n)*A',n,n,n),[3 2 1]);
% 
% temp = reshape(reshape(p3,[],n)*A',n,n,n);
% temp = permute(reshape(reshape(permute(temp,[1 3 2]),[],n)*A',n,n,n),[1 3 2]);
% p3_monomial = permute(reshape(reshape(permute(temp,[3 2 1]),[],n)*A',n,n,n),[3 2 1]);

% The Cayley function
f_cayley = @(s1,s2,t1,t2,z) (cheb1(s1,s2,z).*cheb2(t1,s2,z).*cheb3(t1,t2,z) + cheb2(s1,s2,z).*cheb3(t1,s2,z).*cheb1(t1,t2,z) + cheb3(s1,s2,z).*cheb1(t1,s2,z).*cheb2(t1,t2,z) ...
    - cheb3(s1,s2,z).*cheb2(t1,s2,z).*cheb1(t1,t2,z) - cheb2(s1,s2,z).*cheb1(t1,s2,z).*cheb3(t1,t2,z) - cheb1(s1,s2,z).*cheb3(t1,s2,z).*cheb2(t1,t2,z)) ./ ((s1-t1).*(s2-t2));

% % s1,s2,z,t1,t2
% M = zeros(n,n,n,n,n,3,3);
% 
% M(:,:,:,1,1,1,1) = p1;
% M(:,:,:,1,1,1,2) = p2;
% M(:,:,:,1,1,1,3) = p3;
% 
% M(1,:,:,1,:,2,1) = permute(p1, [2,3,1]);
% M(1,:,:,1,:,2,2) = permute(p2, [2,3,1]);
% M(1,:,:,1,:,2,3) = permute(p3, [2,3,1]);
% 
% M(1,1,:,:,:,3,1) = permute(p1, [3,1,2]);
% M(1,1,:,:,:,3,2) = permute(p2, [3,1,2]);
% M(1,1,:,:,:,3,3) = permute(p3, [3,1,2]);

% % In monomial basis
% % s1,s2,t1,t2,z
% M = zeros(n,n,n,n,n,3,3);
% 
% M(:,:,1,1,:,1,1) = p1_monomial;
% M(:,:,1,1,:,1,2) = p2_monomial;
% M(:,:,1,1,:,1,3) = p3_monomial;
% 
% M(1,:,1,:,:,2,1) = permute(p1_monomial, [2,1,3]);
% M(1,:,1,:,:,2,2) = permute(p2_monomial, [2,1,3]);
% M(1,:,1,:,:,2,3) = permute(p3_monomial, [2,1,3]);
% 
% M(1,1,:,:,:,3,1) = p1_monomial;
% M(1,1,:,:,:,3,2) = p2_monomial;
% M(1,1,:,:,:,3,3) = p3_monomial;

% tic
% cayley_product = five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,1),M(:,:,:,:,:,2,2)),M(:,:,:,:,:,3,3)) + five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,2),M(:,:,:,:,:,2,3)),M(:,:,:,:,:,3,1)) + ...
% five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,3),M(:,:,:,:,:,2,1)),M(:,:,:,:,:,3,2)) - five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,3),M(:,:,:,:,:,2,2)),M(:,:,:,:,:,3,1)) - ...
% five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,2),M(:,:,:,:,:,2,1)),M(:,:,:,:,:,3,3)) - five_variate_multiply(five_variate_multiply(M(:,:,:,:,:,1,1),M(:,:,:,:,:,2,3)),M(:,:,:,:,:,3,2));
% 
% A2 = cayley_divide(cayley_product);
% norm(A2(:))
% A2 = A2(1:2*n-1,1:2*n-1,1:2*n-1,1:2*n-1,:);
% norm(A2(:))
% toc


D1_adj = @(s1,s2,t1,t2,z) [0 -(Dx_cheb2(s1,s2,z).*cheb3(t1,t2,z) - Dx_cheb3(s1,s2,z).*cheb2(t1,t2,z)) (Dx_cheb2(s1,s2,z).*cheb3(t1,s2,z) - Dx_cheb3(s1,s2,z).*cheb2(t1,s2,z)); ...
    0 (Dx_cheb1(s1,s2,z).*cheb3(t1,t2,z) - Dx_cheb3(s1,s2,z).*cheb1(t1,t2,z)) -(Dx_cheb1(s1,s2,z).*cheb3(t1,s2,z) - Dx_cheb3(s1,s2,z).*cheb1(t1,s2,z)); ... 
    0 -(Dx_cheb1(s1,s2,z).*cheb2(t1,t2,z) - Dx_cheb2(s1,s2,z).*cheb1(t1,t2,z)) (Dx_cheb1(s1,s2,z).*cheb2(t1,s2,z) - Dx_cheb2(s1,s2,z).*cheb1(t1,s2,z))];

% D1_M = zeros(n,n,n,n,n,3,3);
% D2_M = zeros(n,n,n,n,n,3,3);
% D12_M = zeros(n,n,n,n,n,3,3);
% 
% D1_M(1:n-1,:,:,:,:,:,:) = M(2:n,:,:,:,:,:,:).*repmat((1:n-1)',1,n,n,n,n,3,3);
% D2_M(:,1:n-1,:,:,:,:,:) = M(:,2:n,:,:,:,:,:).*repmat((1:n-1),n,1,n,n,n,3,3);
% D12_M(1:n-1,1:n-1,:,:,:,:,:) = M(2:n,2:n,:,:,:,:,:).*repmat((1:n-1)',1,n-1,n,n,n,3,3).*repmat((1:n-1),n-1,1,n,n,n,3,3);

% We know the values of the functions at the chebyshev points, so you could make this faster by replacing the evaluation of chebfuns with the already known values 
D_f_cayley = @(s1,s2,t1,t2,z,D_M) trace(adjoint([cheb1(s1,s2,z) cheb2(s1,s2,z) cheb3(s1,s2,z); cheb1(t1,s2,z) cheb2(t1,s2,z) cheb3(t1,s2,z); cheb1(t1,t2,z) cheb2(t1,t2,z) cheb3(t1,t2,z)]) ... 
   * D_M);

D12_f_cayley = @(s1,s2,t1,t2,z,D2_M,D12_M) trace(adjoint([cheb1(s1,s2,z) cheb2(s1,s2,z) cheb3(s1,s2,z); cheb1(t1,s2,z) cheb2(t1,s2,z) cheb3(t1,s2,z); cheb1(t1,t2,z) cheb2(t1,t2,z) cheb3(t1,t2,z)]) ... 
    * D12_M + D1_adj(s1,s2,t1,t2,z) * D2_M);

% Take f_cayley as a 5D function and interpolate it assuming we know the value of the
% function at the Cheby points, not considering how to get these values
% from the nodal values of the functions f1, f2 and f3.

% Where does this number come from? If two highest coeffs are 0, then
% determinant is of degree 0 (in monomial basis implementation at least)
n=2*n-1;
% 5D Chebyshev points
[a,b,c,d,e] = ndgrid(cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)');

% Evaluate the function at the given points (a,b,c,d,e)
f = f_cayley(a(:), b(:), c(:), d(:), e(:));
% Express the values in a grid instead of a list
f = reshape(f,n,n,n,n,n);
toc

tic
D1_vals = zeros(n,n,n,n,n,3,3);
D2_vals = zeros(n,n,n,n,n,3,3);
D12_vals = zeros(n,n,n,n,n,3,3);
% Use chebfun to find the values of the derivative?
for j=1:3
    f_1 = Dx_cheb(j);
    f_2 = Dy_cheb(j);
    f_3 = Dxy_cheb(j);

    temp1 = f_1(a(:), b(:), e(:));
    temp2 = f_2(a(:), b(:), e(:));
    temp3 = f_3(a(:), b(:), e(:));

    % Express the values in a grid instead of a list
    temp1 = reshape(temp1,n,n,n,n,n);
    temp2 = reshape(temp2,n,n,n,n,n);
    temp3 = reshape(temp3,n,n,n,n,n);
    
    D1_vals(:,:,:,:,:,1,j) = temp1;
    D2_vals(:,:,:,:,:,1,j) = temp2;
    D12_vals(:,:,:,:,:,1,j) = temp3;
    
    
    
    temp1 = f_1(c(:), b(:), e(:));
    temp2 = f_2(c(:), b(:), e(:));
    temp3 = f_3(c(:), b(:), e(:));

    % Express the values in a grid instead of a list
    temp1 = reshape(temp1,n,n,n,n,n);
    temp2 = reshape(temp2,n,n,n,n,n);
    temp3 = reshape(temp3,n,n,n,n,n);
    
    D1_vals(:,:,:,:,:,2,j) = temp1;
    D2_vals(:,:,:,:,:,2,j) = temp2;
    D12_vals(:,:,:,:,:,2,j) = temp3;
    
    
    
    temp1 = f_1(c(:), d(:), e(:));
    temp2 = f_2(c(:), d(:), e(:));
    temp3 = f_3(c(:), d(:), e(:));

    % Express the values in a grid instead of a list
    temp1 = reshape(temp1,n,n,n,n,n);
    temp2 = reshape(temp2,n,n,n,n,n);
    temp3 = reshape(temp3,n,n,n,n,n);
    
    D1_vals(:,:,:,:,:,3,j) = temp1;
    D2_vals(:,:,:,:,:,3,j) = temp2;
    D12_vals(:,:,:,:,:,3,j) = temp3;
    
    
    
end
% for i=1:3
%     for j=1:3
%         temp1 = poly_eval_5_variate(D1_M(:,:,:,:,:,i,j),[a(:) b(:) c(:) d(:) e(:)]);
%         temp2 = poly_eval_5_variate(D2_M(:,:,:,:,:,i,j),[a(:) b(:) c(:) d(:) e(:)]);
%         temp3 = poly_eval_5_variate(D12_M(:,:,:,:,:,i,j),[a(:) b(:) c(:) d(:) e(:)]);
%         % Express the values in a grid instead of a list
%         D1_vals(:,:,:,:,:,i,j) = reshape(temp1,n,n,n,n,n);
%         D2_vals(:,:,:,:,:,i,j) = reshape(temp2,n,n,n,n,n);
%         D12_vals(:,:,:,:,:,i,j) = reshape(temp3,n,n,n,n,n);
%     end
% end
toc 
tic
p = cos((2*(1:n)-1)/(2*n)*pi)';
for s1=1:n
   for s2=1:n
       for z=1:n
           for t=1:n
               f(s1,s2,s1,t,z) = D_f_cayley(p(s1),p(s2),p(s1),p(t),p(z),reshape(D1_vals(s1,s2,s1,t,z,:,:),3,3))/(p(s2)-p(t));
               f(s1,s2,t,s2,z) = D_f_cayley(p(s1),p(s2),p(t),p(s2),p(z),reshape(D2_vals(s1,s2,t,s2,z,:,:),3,3))/(p(s1)-p(t));
           end
           f(s1,s2,s1,s2,z) = D12_f_cayley(p(s1),p(s2),p(s1),p(s2),p(z),reshape(D2_vals(s1,s2,s1,s2,z,:,:),3,3),reshape(D12_vals(s1,s2,s1,s2,z,:,:),3,3));
       end
   end
end
toc

% f = zeros(2,2,2,2,2);
% f(:,:,:,:,1) = p(1);
% f(:,:,:,:,2) = p(2);

% Degree is kn-1 for all tau_k
A = cheby_5D_interpolate(f,n);

% In ideal case
% N = factorial(d-1)*n^(d-1); R = reshape(A, N, N);

R = reshape(A, n^2, n^2, n);
% R = reshape(A2, n^2, n^2, []);

z_length = n;
for i = flip(1:n)
    if (norm(R(:,:,i),'fro') < 1e-10); z_length = i-1; else; break; end
end

if (z_length <= 2)
    [V,D] = eig(R(1,:,:),R(2,:,:));
    roots = diag(D)
else
    n = n^2;

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

    C1 = -2*eye(n*(z_length-1));
    C1(end-n+1:end,end-n+1:end) = 2*R(:,:,z_length);
    C1(1:n,1:n) = -eye(n);

    C2 = zeros(size(C1));
    C2(1:end-n,n+1:end) = eye(n*(z_length-2));
    C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)+eye(n*(z_length-2));
    C2(end-n+1:end,end-2*n+1:end-n) = -R(:,:,z_length);
    % C2(end-n+1:end, end-2*n+1:end-n) = -2*eye(n,n);

    % Compute the first rows of the coefficient matrix C2
    D=[];
    for i = 1:z_length-1
        D = [D R(:,:,i)];
    end
    C2(end-n+1:end,:) = C2(end-n+1:end,:)+D;

    % Solve the eigenproblem
    [V,D] = eig(C2,-C1);

    roots = diag(D)
end

% R = reshape(A2, n, n, []);
% 
% % Compute the linearization matrices C1 and C2
% D1 = -eye(n*(z_length-1));
% D1(end-n+1:end,end-n+1:end) = R(:,:,z_length)+rand(n)*1e-50; %flip(flip(A*R(:,:,end)*A',1),2);
% 
% D2 = zeros(size(D1));
% D2(1:end-n,n+1:end) = eye(n*(z_length-2));
% 
% % Compute the first rows of the coefficient matrix C2
% D=[];
% for i = 1:z_length-1
%     D = [D R(:,:,i)];
% end
% D2(end-n+1:end,:) = D;
% 
% % Solve the eigenproblem
% [V,D] = eig(D2,-D1);
% 
% roots = diag(D)
