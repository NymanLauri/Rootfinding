function roots = subdivide_rootfinder(f1,f2,f3)

%How many times subdivision needs to be performed in each variable
%direction
x_divide = 0;
y_divide = 0;
z_divide = 0;

% addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')

input_vals.nodes = [-1,1];

result = subdivide(f1,f2,f3,input_vals,1);
nodes = [result.nodes 1];
% 
% affine = @(x) 1/2*(x+1)*(nodes(2) - nodes(1)) + nodes(1); 
% affine = @(x) 1/2*x*(nodes(2) - nodes(1)) + 1/2*(nodes(2) + nodes(1)); 
% affine = @(x) (x - p1)*nodes(2)/(p2-p1) + (x - p2)*nodes(1)/(p1-p2);

affine = @(x) (x - nodes(1))/(nodes(2) - nodes(1)) - (x - nodes(2))/(nodes(1)-nodes(2));
affine_inv = @(x) 1/2*x*(nodes(2) - nodes(1)) + 1/2*(nodes(2) + nodes(1)); 

f1_local = @(x,y,z) f1(affine(x),y,z);
f2_local = @(x,y,z) f2(affine(x),y,z);
f3_local = @(x,y,z) f3(affine(x),y,z);

result.degrees(1)

roots = trivariate_rootfinder(f1_local,f2_local,f3_local,result.degrees(1)-1);

% roots = affine_inv(roots)

roots;

end