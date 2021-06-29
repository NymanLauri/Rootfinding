clear all
addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')

f = @(x,y,z) sin(0.8*x).*cos(z.*y);
g = @(x,y,z) y.^4;
h = @(x,y,z) z.^4;

input_vals.nodes = [-1,1];
iter_range = 0:10
for iter_max = iter_range
    result = helper2(f,g,h,input_vals,iter_max);
%     result = helper1(f,g,h,input_vals,1,iter_max);
    % subdivide_rootfinder(f,g,h);

    degrees(iter_max+1) = result.degrees(1);
    iter_max
end
% a = toc;

close all
plot(iter_range, degrees); hold on
xlabel('iteration')
ylabel('degree')

fun = @(k,d,d1) sum((d-d1.*k.^iter_range).^2);

% Degree reduction coefficient. Should be less than 0.79. Two ways of
% computing.
x1 = fminbnd(@(k) fun(k,degrees,degrees(1)), 0,1)
x2 = (degrees(end)/degrees(1))^(1/iter_range(end))

plot(iter_range, degrees(1)*x1.^iter_range)


function grid_out = helper2(f1,f2,f3,grid_in,iter_max)
interval = [-1 -1+2*0.5^iter_max];
f1_cheb = chebfun3(@(x,y,z) f1(x,y,z), [interval -1 1 -1 1]);
[n f1_2 f1_3] = length(f1_cheb);
grid_out.degrees = [n];

end


function grid_out = helper1(f1,f2,f3,grid_in,iteration, iter_max)
% Max amount of subdivisions is 10 for now
if  iteration <= iter_max
    input_first.nodes = [grid_in.nodes(1) 1/2*(grid_in.nodes(1)+grid_in.nodes(2))];
    grid_first = helper1(f1,f2,f3,input_first,iteration + 1, iter_max);
    
    input_last.nodes = [1/2*(grid_in.nodes(1)+grid_in.nodes(2)) grid_in.nodes(2)];
    grid_last = helper1(f1,f2,f3,input_last,iteration + 1, iter_max);
    
    grid_out.nodes = [grid_first.nodes grid_last.nodes];
    grid_out.degrees = [grid_first.degrees grid_last.degrees];
else
    f1_cheb = chebfun3(@(x,y,z) f1(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);
    f2_cheb = chebfun3(@(x,y,z) f2(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);
    f3_cheb = chebfun3(@(x,y,z) f3(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1]);

    [f1_1 f1_2 f1_3] = length(f1_cheb);
    [f2_1 f2_2 f2_3] = length(f2_cheb);
    [f3_1 f3_2 f3_3] = length(f3_cheb);
    n = max([f1_1, f2_1, f3_1]);  
    
    grid_out.nodes = grid_in.nodes(1);
    grid_out.degrees = [n];
end

end
