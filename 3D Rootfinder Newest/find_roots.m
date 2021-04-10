function roots = find_roots(f1,f2,f3,n)

z_roots = trivariate_rootfinder(f1,f2,f3,n);
roots = [];

tic

% TODO: Choose in which order you solve z, y and x based on the degrees to
% minimize running time
for i=1:size(z_roots)
    if (isinf(z_roots(i)) || ~isreal(z_roots(i)) ); continue; end
    h1 = @(x,y) f1(x,y,z_roots(i));
    h2 = @(x,y) f2(x,y,z_roots(i));
%     h3 = @(x,y) f3(x,y,z_roots(i));
    y_roots = bivariate_rootfinder(h1,h2,n);

    for j=1:size(y_roots)
        if (isinf(y_roots(j)) || ~isreal(y_roots(j))); continue; end
        g1 = @(x) h1(x,y_roots(j));
%         g2 = @(x) h2(x,y_roots(j));
%         g3 = @(x) h3(x,y_roots(j));
        y_roots;
        x_roots = univariate_rootfinder(g1,n);
        treshold = 1e1;
        sols = (abs(f1(x_roots,y_roots(j),z_roots(i))) < treshold) & (abs(f2(x_roots,y_roots(j),z_roots(i))) < treshold) & (abs(f3(x_roots,y_roots(j),z_roots(i))) < treshold);
        k = find(sols);
        if (size(k,1) ~= 0)
            roots = [roots; x_roots(k) repmat(y_roots(j),size(k,1),1) repmat(z_roots(i),size(k,1),1)];
        end
    end
end

roots;
toc

treshold=1e-12;
syms x y z;
oldroots=roots;
J=matlabFunction([diff(f1,x) diff(f2,x) diff(f3,x);diff(f1,y) diff(f2,y) diff(f3,y);diff(f1,z) diff(f2,z) diff(f3,z)]);
for ns=1:size(roots,1)
    g=roots(ns,:);
        for steps=1:100
            g=g-[f1(g(1),g(2),g(3)) f2(g(1),g(2),g(3)) f3(g(1),g(2),g(3))]/J(g(1),g(2),g(3));
            if (abs(f1(g(1),g(2),g(3))) < treshold) & (abs(f2(g(1),g(2),g(3))) < treshold) & (abs(f3(g(1),g(2),g(3))) < treshold)
                break;
            end
        end
    roots(ns,:)=g;
end

% roots;
% polishment=abs(oldroots-roots)

% treshold=1e-5;
% sols = (abs(f1(roots(:,1),roots(:,2),roots(:,3))) < treshold) & (abs(f2(roots(:,1),roots(:,2),roots(:,3))) < treshold) & (abs(f3(roots(:,1),roots(:,2),roots(:,3))) < treshold);
% k = find(sols);
% roots_final = roots(k,:);

end