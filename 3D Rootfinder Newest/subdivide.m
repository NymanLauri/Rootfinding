function grid_out = subdivide(f1,f2,f3,grid_in,iteration)

addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
f1_cheb = chebfun3(@(x,y,z) f1(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);
f2_cheb = chebfun3(@(x,y,z) f2(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);
f3_cheb = chebfun3(@(x,y,z) f3(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);

[f1_1 f1_2 f1_3] = length(f1_cheb);
[f2_1 f2_2 f2_3] = length(f2_cheb);
[f3_1 f3_2 f3_3] = length(f3_cheb);
n = max([f1_1, f2_1, f3_1]);

% Max amount of subdivisions is 10 for now
if n > 5 && iteration < 10
    input_first.nodes = [grid_in.nodes(1) 1/2*(grid_in.nodes(1)+grid_in.nodes(2))];
    grid_first = subdivide(f1,f2,f3,input_first,iteration + 1);
    
    input_last.nodes = [1/2*(grid_in.nodes(1)+grid_in.nodes(2)) grid_in.nodes(2)];
    grid_last = subdivide(f1,f2,f3,input_last,iteration + 1);
    
    grid_out.nodes = [grid_first.nodes grid_last.nodes];
    grid_out.degrees = [grid_first.degrees grid_last.degrees];
else
    grid_out.nodes = grid_in.nodes(1);
    grid_out.degrees = [n];
end

end