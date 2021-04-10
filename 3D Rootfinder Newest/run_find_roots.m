% % "Bad" example
% u = 1;
% n=5;
% q1 = @(x,y,z) x.^n + u * (x+y+z)*sqrt(1/3) + z.^n + y.^n;
% q2 = @(x,y,z) y.^n + u * (x*sqrt(2/3) - y*sqrt(1/6) - z*sqrt(1/6)) + z.^n + x.^n;
% q3 = @(x,y,z) z.^n + u * (y - z)*sqrt(1/2);
% % n = n+1;
% 
% x0 = 0.1;
% y0 = 0.2;
% z0 = 0.3;
% f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
% f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
% f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);
% 
% % n = 5;
% % f1 = @(x,y,z) (y-1/sqrt(2)).^n.*(1+z).^2 + x.^2 + y.^2 - 1;
% % f2 = @(x,y,z) x-y + z.^5-0.5^5;
% % f3 = @(x,y,z) z.^n - 0.5^n + x - 1/sqrt(2);

% Normal example
% n = 2;
% q1 = @(x,y,z) ((x-0.5).^2 + (y-0.5).^2 + z.^2 - 0.5);
% q2 = @(x,y,z) ((x+0.5).^2 + (y-0.5).^2 + z.^2 - 0.5);
% q3 = @(x,y,z) (x.^2 + y.^2 + z.^2 - 0.5^2);
% 
% x0 = 0.1; %The other root is more inaccurate, and changing to +0.1 makes it even more so. Why?
% y0 = 0.0;
% z0 = 0.0;
% f1 = @(x,y,z) q1(x-x0,y-y0,z-z0);
% f2 = @(x,y,z) q2(x-x0,y-y0,z-z0);
% f3 = @(x,y,z) q3(x-x0,y-y0,z-z0);
% 
% % roots - [x0, 0.25, -sqrt(3)/4; x0, 0.25, sqrt(3)/4]

% Non-polynomial example
n=6;
f1 = @(x,y,z) cos(2*pi*x).*cos(2*pi*y).*cos(2*pi*z);
f2 = @(x,y,z) y; %sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
f3 = @(x,y,z) x.^2 + y.^2 + z.^2 - 1;
% sqrt(15)/4, 0.25, 0.75, then sin(arccos(0.75))=sqrt(7)/4
% vals = [sqrt(15)/4, 0.25, 0.75, sqrt(7)/4] % Each has plus and minus
exact = [0.25 0 sqrt(15)/4; 0.25 0 -sqrt(15)/4; -0.25 0 sqrt(15)/4; -0.25 0 -sqrt(15)/4; ...
    0.75 0 sqrt(7)/4; 0.75 0 -sqrt(7)/4; -0.75 0 sqrt(7)/4; -0.75 0 -sqrt(7)/4; ... 
    sqrt(15)/4 0 0.25; sqrt(15)/4 0 -0.25; -sqrt(15)/4 0 0.25; -sqrt(15)/4 0 -0.25; ...
    sqrt(7)/4 0 0.75; sqrt(7)/4 0 -0.75; -sqrt(7)/4 0 0.75; -sqrt(7)/4 0 -0.75];
% min(abs(roots(:,3) - vals(1))) % largest error of order of magnitude 1e-07 with
% not interpolating (n=7) as well as with interpolating.

roots = find_roots(f1,f2,f3,n);

% syms x y z;
% oldroots=roots;
% J=matlabFunction([diff(f1,x) diff(f2,x) diff(f3,x);diff(f1,y) diff(f2,y) diff(f3,y);diff(f1,z) diff(f2,z) diff(f3,z)]);
% for ns=1:size(roots,1)
%     g=roots(ns,:);
%         for steps=1:10
%             g=g-[f1(g(1),g(2),g(3)) f2(g(1),g(2),g(3)) f3(g(1),g(2),g(3))]/J(g(1),g(2),g(3));
%         end
%     roots(ns,:)=g;
% end
% roots;
% polishment=abs(oldroots-roots)
% 
% treshold=1e-5;
% sols = (abs(f1(roots(:,1),roots(:,2),roots(:,3))) < treshold) & (abs(f2(roots(:,1),roots(:,2),roots(:,3))) < treshold) & (abs(f3(roots(:,1),roots(:,2),roots(:,3))) < treshold);
% k = find(sols);
% roots_final = roots(k,:);

for i=1:size(exact,1)
    min(vecnorm((roots - exact(i,:))'))
end

