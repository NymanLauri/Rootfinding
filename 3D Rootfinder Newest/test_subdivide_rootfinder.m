clear all
addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
% addpath('C:\Users\Lauri\Documents\MATLAB\Rootfinding\chebfun-master')

% Root at z = y = 1/sqrt(2) and when sine=1
f = @(x,y,z) 1/4*sin(2*pi*2.1*x) + y.^4;
g = @(x,y,z) z.^2 + y.^2 - 1;
h = @(x,y,z) z-y;

subdivide_rootfinder(f,g,h);