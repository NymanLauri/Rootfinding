function roots = subdivide_rootfinder(f1,f2,f3)

%How many times subdivision needs to be performed in each variable
%direction
x_divide = 0;
y_divide = 0;
z_divide = 0;

% addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')

input_vals.nodes = [-1,1];

result = subdivide(f1,f2,f3,input_vals,1);
nodes = [result.nodes 1]
% 
% affine = @(x) 1/2*(x+1)*(nodes(2) - nodes(1)) + nodes(1); 
% affine = @(x) 1/2*x*(nodes(2) - nodes(1)) + 1/2*(nodes(2) + nodes(1)); 
% affine = @(x) (x - p1)*nodes(2)/(p2-p1) + (x - p2)*nodes(1)/(p1-p2);

affine = @(x,a,b) 1/2*x*(b - a) + 1/2*(b + a); 
% affine_inv = @(x,a,b) (x - a)/(b - a) - (x - b)/(a-b);

sum(result.roots)
roots = [];
for i=1:length(result.degrees)
    if result.roots(i)
        result.degrees(i)
        f1_local = @(x,y,z) f1(affine(x,nodes(i),nodes(i+1)),y,z);
        f2_local = @(x,y,z) f2(affine(x,nodes(i),nodes(i+1)),y,z);
        f3_local = @(x,y,z) f3(affine(x,nodes(i),nodes(i+1)),y,z);
%         roots = trivariate_rootfinder(f1_local,f2_local,f3_local,result.degrees(i));
        roots_local = find_roots(f1_local,f2_local,f3_local,result.degrees(i));
%         roots_local = find_roots(f1,f2,f3,result.degrees(i));
        if ~isempty(roots_local); roots_local = [affine(roots_local(:,1),nodes(i),nodes(i+1)) roots_local(:,2:3)]; end;
        roots = [roots; roots_local];
    end
end
% roots = affine_inv(roots)

roots;

end