clear all
addpath('..\chebfun-master')

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

% Solving the eigenvalue problem happens takes less time than the 
% evaluation and interpolation of the Cayley function since the degrees 
% z are low. 
roots = subdivide_rootfinder(f,g,h);
[f(roots(:,1),roots(:,2),roots(:,3)) g(roots(:,1),roots(:,2),roots(:,3)) h(roots(:,1),roots(:,2),roots(:,3))];


% Without subdivision for running time comparison, becomes infeaseible
precision = 1e-12;

f_cheb = chebfun3t(@(x,y,z) f(x,y,z), 'eps', precision);
g_cheb = chebfun3t(@(x,y,z) g(x,y,z), 'eps', precision);
h_cheb = chebfun3t(@(x,y,z) h(x,y,z), 'eps', precision);

[f_1 f_2 f_3] = length(f_cheb);
[g_1 g_2 g_3] = length(g_cheb);
[h_1 h_2 h_3] = length(h_cheb);
n = max([f_1, f_2, f_3, g_1, g_2, g_3, h_1, h_2, h_3])-1

roots2 = find_roots(f,g,h,n);
[f(roots2(:,1),roots2(:,2),roots2(:,3)) g(roots2(:,1),roots2(:,2),roots2(:,3)) h(roots2(:,1),roots2(:,2),roots2(:,3))];


