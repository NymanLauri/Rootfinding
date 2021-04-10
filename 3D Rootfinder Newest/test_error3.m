% % Inaccuracy of 0.1
% f1 = @(x,y,z) y.^2 + x.^2 - 1;
% f2 = @(x,y,z) (-z+0.5)+x-1/sqrt(2);
% f3 = @(x,y,z) (z-0.5).*(2*x.^2 + y.^2 - 1); %Removing coefficient 2 creates a singular jacobian and instability
% 
% n = 2;
% 
% % Root should be 0.5
% roots = trivariate_rootfinder(f1,f2,f3,n)



% % Random tests that can be removed
f1 = @(x,y,z) (y-0.7).^2 + z.*x.^2 - 1;
f2 = @(x,y,z) (y+0.5).^2 + (x-1.5).^2 - 1;
f3 = @(x,y,z) y.*(y+0.6).^2 + (x+1).^2 - 1;

% f1 = @(x,y,z) (y-0.7).^2 + z.*x.^2 - 1;
% f2 = @(x,y,z) 0;
% f3 = @(x,y,z) 0;

n = 3;

roots = trivariate_rootfinder(f1,f2,f3,n)



