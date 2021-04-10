% Error behaviour of different "bad systems"
vals1 = [];
vals2 = [];
vals3 = [];
for u = flip(10.^(-16:-1))
    % Orthogonal matrix 1
    q1 = @(x,y,z) x.^2 + u * (x+y+z)*sqrt(1/3);
    q2 = @(x,y,z) y.^2 + u * (x*sqrt(2/3) - y*sqrt(1/6) - z*sqrt(1/6));
    q3 = @(x,y,z) z.^2 + u * (y - z)*sqrt(1/2);
    n=2;
        
    x0 = 2*rand-1;
    y0 = 2*rand-1;
    z0 = 2*rand-1;

    f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
    f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
    f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);

    roots = trivariate_rootfinder(f1,f2,f3,n);
    vals1 = [vals1 min(abs(roots - z0))];
    
    
    % Orthogonal matrix 2
    q1 = @(x,y,z) x.^2 + u*(x);
    q2 = @(x,y,z) y.^2 + u*(sqrt(1/2)*y + sqrt(1/2)*z);
    q3 = @(x,y,z) z.^2 + u*(sqrt(1/2)*y - sqrt(1/2)*z);
    n=2;
    
    f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
    f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
    f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);
    
    roots = trivariate_rootfinder(f1,f2,f3,n);
    vals2 = [vals2 min(abs(roots - z0))];
    
    % Orthogonal matrix 3
    q1 = @(x,y,z) x.^2 + u*(sqrt(3)/2*x + sqrt(3)/4*y + 1/4*z);
    q2 = @(x,y,z) y.^2 + u*(-1/2*x + 3/4*y + sqrt(3)/4*z);
    q3 = @(x,y,z) z.^2 + u*(-1/2*y + sqrt(3)/2*z);
    n=2;
    
    f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
    f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
    f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);
    
    roots = trivariate_rootfinder(f1,f2,f3,n);
    vals3 = [vals3 min(abs(roots - z0))];
end
close all
figure
subplot(3,1,1)
loglog(flip(10.^(-16:-1)),vals1, 'x')
xlabel('u')
ylabel('inaccuracy of the root')
title('First')
subplot(3,1,2)
loglog(flip(10.^(-16:-1)),vals2, 'x')
xlabel('u')
ylabel('inaccuracy of the root')
title('Second')
subplot(3,1,3)
loglog(flip(10.^(-16:-1)),vals3, 'x')
xlabel('u')
ylabel('inaccuracy of the root')
title('Third')

figure
loglog(flip(10.^(-16:-1)),vals1, '-x',flip(10.^(-16:-1)),vals2, '-x',flip(10.^(-16:-1)),vals3, '-x')
legend('First','Second','Third')
xlabel('u')
ylabel('inaccuracy of the root')


% J = [sqrt(1/3)*u sqrt(1/3)*u sqrt(1/3)*u; sqrt(2/3)*u -sqrt(1/6)*u -sqrt(1/6)*u; 0 sqrt(1/2)*u -sqrt(1/2)*u];
% norm(inv(J));
