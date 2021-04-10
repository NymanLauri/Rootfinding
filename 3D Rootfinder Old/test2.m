% fun = @(a,b,c,d,e) a.^2;
% 
% n=3
% 
% [a,b,c,e,g] = ndgrid(cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)',cos((2*(1:n)-1)/(2*n)*pi)');
% 
% % Evaluate the function at the given points (a,b,c,d,e)
% f = fun(a(:), b(:), c(:), e(:), g(:));
% % Express the values in a grid instead of a list
% f = reshape(f,n,n,n,n,n);
% 
% cheby_5D_interpolate(f,n)


% n = 3;
% R = rand(n,n,n,n,n);
% R = reshape(R, n^2, n^2, n);
% R(:,:,end) = 0;
% n = n^2;

n = 3;
R = zeros(n,n,n,n,n);
% R(1:2,1:2,1:2,1:2,1:2) = rand(2,2,2,2,2);
R = reshape(R, n^2, n^2, n);
n = n^2;

R([1 2],[1 2],1) = rand(2);
% R([1 2 4 5],[1 2 4 5],2) = -rand(4);
R(1:2,1:2,2) = eye(2);

% R = rand(n,n,sqrt(n))
% R(:,:,end)=0;

% R([1 4],[1 3],1) = rand(2);
% R([1 3],[1 3],2) = -rand(2);
% R(:,:,2) = -eye(n);

A = cheby_basis(n);

% % Compute the linearization matrices C1 and C2
% C1 = -2*eye(n*(size(R,3)-1));
% C1(end-n+1:end,end-n+1:end) = 2*R(:,:,end);
% C1(1:n,1:n) = -eye(n);
% 
% C2 = zeros(size(C1));
% C2(1:end-n,n+1:end) = eye(n*(size(R,3)-2));
% C2(n+1:end,1:end-n) = C2(n+1:end,1:end-n)+eye(n*(size(R,3)-2));
% C2(end-n+1:end,end-2*n+1:end-n) = -R(:,:,end);
% % C2(end-n+1:end, end-2*n+1:end-n) = -2*eye(n,n);
% 
% % Compute the first rows of the coefficient matrix C2
% D=[];
% for i = 1:size(R,3)-1
%     D = [D R(:,:,i)];
% end
% C2(end-n+1:end,:) = C2(end-n+1:end,:)+D;
% % C2 = C2;
% 
% % Compute the linearization matrices C1 and C2
% D1 = eye(n*(size(R,3)-1));
% D1(end-n+1:end,end-n+1:end) = R(:,:,end); %+1e-10*eye(size(R(:,:,end))); %flip(flip(A*R(:,:,end)*A',1),2);
% 
% D2 = zeros(size(D1));
% D2(1:end-n,n+1:end) = eye(n*(size(R,3)-2));
% 
% % Compute the first rows of the coefficient matrix C2
% D=[];
% for i = 1:size(R,3)-1
%     D = [D R(:,:,i)];
% end
% D2(end-n+1:end,:) = -D;

% Compute the linearization matrices C1 and C2
D1 = eye(n*(size(R,3)-1));
D1(1:n,1:n) = R(:,:,end); %+1e-10*eye(size(R(:,:,end))); %flip(flip(A*R(:,:,end)*A',1),2);

D2 = zeros(size(D1));
D2(n+1:end,1:end-n) = eye(n*(size(R,3)-2));

% Compute the first rows of the coefficient matrix C2
D=[];
for i = flip(1:size(R,3)-1)
    D = [D R(:,:,i)];
end
D2(1:n,:) = -D;

last = n*(size(R,3)-1);
idx = 1:last;
idx([1 2 3 4 5]) = idx([1 2 4 5 3]);

idx2 = idx;
idx2([1 2 3 last-1 last]) = idx([last 2 last-1 1 3]);


% eig(C2(idx2,idx2),-C1(idx2,idx2))
% eig(C2,C1)
eig(D2,D1+1e-250*rand(size(D1)))
% D2 = D2(:,idx2);
% eig(D2+1e-250*rand(size(D1)),D1(:,idx2)+1e-250*rand(size(D1)))
eig(R(:,:,1),-R(:,:,2))
