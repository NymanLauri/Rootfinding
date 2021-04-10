% The effect of the order of first terms
vals1 = [];
vals2 = [];
vals3 = [];
for u = flip(10.^(-16:-1))
    % Original bad example
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
       
    % Original with order 4
    q1 = @(x,y,z) x.^4 + u * (x+y+z)*sqrt(1/3);
    q2 = @(x,y,z) y.^4 + u * (x*sqrt(2/3) - y*sqrt(1/6) - z*sqrt(1/6));
    q3 = @(x,y,z) z.^4 + u * (y - z)*sqrt(1/2);
    n=4;
    
    f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
    f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
    f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);
    
    roots = trivariate_rootfinder(f1,f2,f3,n);
    vals2 = [vals2 min(abs(roots - z0))];
    
    % Original example with order 6
    q1 = @(x,y,z) x.^6 + u * (x+y+z)*sqrt(1/3);
    q2 = @(x,y,z) y.^6 + u * (x*sqrt(2/3) - y*sqrt(1/6) - z*sqrt(1/6));
    q3 = @(x,y,z) z.^6 + u * (y - z)*sqrt(1/2);
    n=6;
    
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
title('x^2')
subplot(3,1,2)
loglog(flip(10.^(-16:-1)),vals2, 'x')
xlabel('u')
ylabel('inaccuracy of the root')
title('x^4')
subplot(3,1,3)
loglog(flip(10.^(-16:-1)),vals3, 'x')
xlabel('u')
ylabel('inaccuracy of the root')
title('x^6')

figure
loglog(flip(10.^(-16:-1)),vals1, '-x',flip(10.^(-16:-1)),vals2, '-x',flip(10.^(-16:-1)),vals3, '-x')
legend('x^2','x^4','x^6')
xlabel('u')
ylabel('inaccuracy of the root')