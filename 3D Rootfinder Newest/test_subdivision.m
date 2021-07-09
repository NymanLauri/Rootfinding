clear all
addpath('C:\Matlab\MATLAB\Rootfinding\chebfun-master')
% addpath('C:\Users\Lauri\Documents\MATLAB\Rootfinding\chebfun-master')

f = @(x,y,z) sin(2*pi*2.1*x).*cos(z.*y);
g = @(x,y,z) cos(2*pi*y);
h = @(x,y,z) cos(2*pi*z);

input_vals.nodes = [-1,1];
iter_range = 0:8
for iter_max = iter_range
%     result = helper2(f,g,h,input_vals,iter_max);
    result = helper1(f,g,h,input_vals,1,iter_max);
%     result = subdivide(f,g,h,input_vals,1);
    % subdivide_rootfinder(f,g,h);

    degrees(iter_max+1) = result.degrees(1);
    iter_max
end
% result = subdivide(f,g,h,input_vals,1)
% a = toc;

% Fraction of subdomains where roots can exist
frac = sum(result.roots)/length(result.roots)

% Desired coeff or degree reduction with this value of frac
exp((-log(frac^(1/iter_range(end))*2))/3)

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
legend('Actual','Fit')

function grid_out = helper2(f1,f2,f3,grid_in,iter_max)
precision = eps;
interval = [-1 -1+2*0.5^iter_max];
f1_cheb = chebfun3(@(x,y,z) f1(x,y,z), [interval -1 1 -1 1], 'eps', precision);
[n f1_2 f1_3] = length(f1_cheb);
grid_out.degrees = [n-1];

end


function grid_out = helper1(f1,f2,f3,grid_in,iteration, iter_max)
precision = 1e-12;
% Max amount of subdivisions is 10 for now
if  iteration <= iter_max
    input_first.nodes = [grid_in.nodes(1) 1/2*(grid_in.nodes(1)+grid_in.nodes(2))];
    grid_first = helper1(f1,f2,f3,input_first,iteration + 1, iter_max);
    
    input_last.nodes = [1/2*(grid_in.nodes(1)+grid_in.nodes(2)) grid_in.nodes(2)];
    grid_last = helper1(f1,f2,f3,input_last,iteration + 1, iter_max);
    
    grid_out.nodes = [grid_first.nodes grid_last.nodes];
    grid_out.degrees = [grid_first.degrees grid_last.degrees];
    grid_out.roots = [grid_first.roots grid_last.roots];
else
    f1_cheb = chebfun3t(@(x,y,z) f1(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);
    f2_cheb = chebfun3t(@(x,y,z) f2(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);
    f3_cheb = chebfun3t(@(x,y,z) f3(x,y,z), [grid_in.nodes(1) grid_in.nodes(2) -1 1 -1 1], 'eps', precision);

    [f1_1 f1_2 f1_3] = length(f1_cheb);
    [f2_1 f2_2 f2_3] = length(f2_cheb);
    [f3_1 f3_2 f3_3] = length(f3_cheb);
    n = max([f1_1, f2_1, f3_1])-1;  
    
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
