function roots = subdivide_rootfinder(f1,f2,f3)

addpath('..\chebfun-master')

input_vals.nodes = [-1,1];

% Subdividing in each variable separately
disp('Time taken for performing the subdivision:')
tic
x_division = subdivide(f1,f2,f3,input_vals,1);
x_nodes = [x_division.nodes 1];

y_division = subdivide(@(x,y,z) f1(z,x,y),@(x,y,z) f2(z,x,y),@(x,y,z) f3(z,x,y),input_vals,1);
y_nodes = [y_division.nodes 1];

z_division = subdivide(@(x,y,z) f1(y,z,x),@(x,y,z) f2(y,z,x),@(x,y,z) f3(y,z,x),input_vals,1);
z_nodes = [z_division.nodes 1];
toc

affine = @(x,a,b) 1/2*x*(b - a) + 1/2*(b + a); 

X = ['Amount of subdomains without additional step: ', num2str(sum(x_division.roots)*sum(y_division.roots)*sum(z_division.roots))];
disp(X)


roots = [];
xyz_roots=zeros(length(x_division.degrees),length(y_division.degrees),length(z_division.degrees));

disp('Time taken for the additional step:')
tic
for i=1:length(x_division.degrees)
    if x_division.roots(i)
        for j=1:length(y_division.degrees)
            if y_division.roots(j)
                for k=1:length(z_division.degrees)
                    if z_division.roots(k)
                        % Checking if solutions can exist in each
                        % subdomain. This should not be the bottle neck.
                        f1_cheb = chebfun3t(@(x,y,z) f1(x,y,z), [x_nodes(i) x_nodes(i+1) y_nodes(j) y_nodes(j+1) z_nodes(k) z_nodes(k+1)]);
                        f2_cheb = chebfun3t(@(x,y,z) f2(x,y,z), [x_nodes(i) x_nodes(i+1) y_nodes(j) y_nodes(j+1) z_nodes(k) z_nodes(k+1)]);
                        f3_cheb = chebfun3t(@(x,y,z) f3(x,y,z), [x_nodes(i) x_nodes(i+1) y_nodes(j) y_nodes(j+1) z_nodes(k) z_nodes(k+1)]);
                        
                        xyz_roots(i,j,k) = (2-eps*10)*abs(f1_cheb.coeffs(1,1,1))<=sum(sum(sum(abs(f1_cheb.coeffs)))) && (2-eps*10)*abs(f2_cheb.coeffs(1,1,1))<=sum(sum(sum(abs(f2_cheb.coeffs)))) ...
                            && (2-eps*10)*abs(f3_cheb.coeffs(1,1,1))<=sum(sum(sum(abs(f3_cheb.coeffs))));
                    end
                end
            end
        end
    end
end
toc


Y = ['Amount of subdomains with additional step: ', num2str(sum(xyz_roots, 'all'))];
disp(Y)

% keyboard 

for i=1:length(x_division.degrees)
    if x_division.roots(i)
        for j=1:length(y_division.degrees)
            if y_division.roots(j)
                for k=1:length(z_division.degrees)
                    if z_division.roots(k)
                        [x_division.degrees(i) y_division.degrees(j) z_division.degrees(k)];
                        
                        if xyz_roots(i,j,k)                       
                            f1_local = @(x,y,z) f1(affine(x,x_nodes(i),x_nodes(i+1)),affine(y,y_nodes(j),y_nodes(j+1)),affine(z,z_nodes(k),z_nodes(k+1)));
                            f2_local = @(x,y,z) f2(affine(x,x_nodes(i),x_nodes(i+1)),affine(y,y_nodes(j),y_nodes(j+1)),affine(z,z_nodes(k),z_nodes(k+1)));
                            f3_local = @(x,y,z) f3(affine(x,x_nodes(i),x_nodes(i+1)),affine(y,y_nodes(j),y_nodes(j+1)),affine(z,z_nodes(k),z_nodes(k+1)));
                            roots_local = find_roots(f1_local,f2_local,f3_local,x_division.degrees(i));
                            if ~isempty(roots_local) 
                                roots_local = [affine(roots_local(:,1),x_nodes(i),x_nodes(i+1)) affine(roots_local(:,2),y_nodes(j),y_nodes(j+1)) affine(roots_local(:,3),z_nodes(k),z_nodes(k+1))]; 
                            end
                            roots = [roots; roots_local];
                        else
%                             [i j k]
                        end
                    end
                end
            end
        end
    end
end


end