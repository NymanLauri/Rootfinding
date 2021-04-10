
% % Only a couple of spurious solutions
% f1 = @(x,y,z) x+y;
% f2 = @(x,y,z) x-y;
% f3 = @(x,y,z) (z-0.123456789);
% 
% n = 3;


% % Spurious solutions
% f1 = @(x,y,z) x+y.*z-1;
% f2 = @(x,y,z) z.*(x-y);
% f3 = @(x,y,z) (z.^3-0.123456789^3);
% 
% n = 3;


% % Spurious solutions
% f1 = @(x,y,z) x+y.*z.^3-1;
% f2 = @(x,y,z) z.*(x-y);
% f3 = @(x,y,z) (z.^4-0.123456789^4);
% 
% n = 4;


% Spurious solutions
f1 = @(x,y,z) y.^2 + x.^2 - 1;
f2 = @(x,y,z) (-z+0.5)+x-1/sqrt(2);
f3 = @(x,y,z) (z-0.5).*(2*x.^2 + y.^2 - 1); %Removing coefficient 2 creates a singular jacobian and instability

n = 6;


roots = trivariate_rootfinder(f1,f2,f3,n)