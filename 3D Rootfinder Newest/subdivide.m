function grid_out = subdivide(f1,f2,f3,grid_in,iteration)
precision = 1e-12;

% addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
f1_cheb = chebfun3t(@(x,y,z) f1(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);
f2_cheb = chebfun3t(@(x,y,z) f2(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);
f3_cheb = chebfun3t(@(x,y,z) f3(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);

[f1_1 f1_2 f1_3] = length(f1_cheb);
[f2_1 f2_2 f2_3] = length(f2_cheb);
[f3_1 f3_2 f3_3] = length(f3_cheb);
n = max([f1_1, f2_1, f3_1])-1;

%Avoid halving the domain exactly in order to not subdivide the domain at a
%root
small_val = 1e-6; 

% Max amount of subdivisions is 10 for now
if n > 5 && iteration <= 10
    input_first.nodes = [grid_in.nodes(1) 1/2*(grid_in.nodes(1)+grid_in.nodes(2))+small_val];
    grid_first = subdivide(f1,f2,f3,input_first,iteration + 1);
    
    input_last.nodes = [1/2*(grid_in.nodes(1)+grid_in.nodes(2))+small_val grid_in.nodes(2)];
    grid_last = subdivide(f1,f2,f3,input_last,iteration + 1);
    
    grid_out.nodes = [grid_first.nodes grid_last.nodes];
    grid_out.degrees = [grid_first.degrees grid_last.degrees];
    grid_out.roots = [grid_first.roots grid_last.roots];
else
    grid_out.nodes = grid_in.nodes(1);
    grid_out.degrees = [n];
    
    if (2-eps*10)*abs(f1_cheb.coeffs(1,1,1))>sum(sum(sum(abs(f1_cheb.coeffs)))) || (2-eps*10)*abs(f2_cheb.coeffs(1,1,1))>sum(sum(sum(abs(f2_cheb.coeffs)))) ...
        || (2-eps*10)*abs(f3_cheb.coeffs(1,1,1))>sum(sum(sum(abs(f3_cheb.coeffs))))
        grid_out.roots = [0];
    else 
        grid_out.roots = [1];
    end
end

end