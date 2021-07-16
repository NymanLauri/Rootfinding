clear all
addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
% addpath('C:\Users\Lauri\Documents\MATLAB\Rootfinding\chebfun-master')

% % Root at z = y = 1/sqrt(2) and when sine=-1
% f = @(x,y,z) 1/4*sin(pi*x) + y.^4;
% g = @(x,y,z) z.^2 + y.^2 - 1;
% h = @(x,y,z) z-y;

% % Root at z = y = 1/sqrt(2) and when sine=-1
% f = @(x,y,z) 1/4*sin(pi*x) + cos(y.^4*pi).^4;
% g = @(x,y,z) z.^2 + y.^2 - 1;
% h = @(x,y,z) z-y;

% Root at z = y = 1/sqrt(2) and when sine=-1
f = @(x,y,z) 1/4*sin(pi*x) - 1/4*(sqrt(2)*y-2).^(-1);
g = @(x,y,z) z.^2 + y.^2 - 1;
h = @(x,y,z) z-y;


roots = subdivide_rootfinder(f,g,h);
[f(roots(:,1),roots(:,2),roots(:,3)) g(roots(:,1),roots(:,2),roots(:,3)) h(roots(:,1),roots(:,2),roots(:,3))];

