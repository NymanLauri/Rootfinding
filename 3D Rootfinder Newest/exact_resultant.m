function ez = exact_resultant(QQ, x0_double, y0_double, z0_double, u_input, U, V)
% clear all
syms x y z u q11 q12 q13 q21 q22 q23 q31 q32 q33 s1 t1 s2 t2 x0 y0 z0
 
%here I compute an algebraic formula for M, the Cayley resultant matpoly
Q=[q11 q12 q13;q21 q22 q23;q31 q32 q33];
vf=[(x-x0)^2;(y-y0)^2;(z-z0)^2]+u*Q*[x-x0;y-y0;z-z0];
f=vf(1);g=vf(2);h=vf(3);
CC=det([[subs(f,{x y},{s1 s2});subs(f,{x y},{t1 s2});subs(f,{x y},{t1 t2})] [subs(g,{x y},{s1 s2});subs(g,{x y},{t1 s2});subs(g,{x y},{t1 t2})] [subs(h,{x y},{s1 s2});subs(h,{x y},{t1 s2});subs(h,{x y},{t1 t2})]]);
CC=simplify(CC/((s1-t1)*(s2-t2)));
one=subs([diff(diff(CC,t1),t2) diff(CC,t1) diff(CC,t2) CC],{t1 t2},{0 0});
M=subs([diff(diff(one,s1),s2);diff(one,s1);diff(one,s2);one],{s1 s2},{0 0});
M=simplify(M);
 
%now I compute the coeffs in Cheby basis and the linearization
I=M^0;O=0*I;
moncoeffs=subs([diff(M,z,2)/2 diff(M,z) M],z,0);
chebcoeffs=moncoeffs*[I/2 O I/2;O I O;O O I];
P2=chebcoeffs(:,1:4);P1=chebcoeffs(:,5:8);P0=chebcoeffs(:,9:12);

P0=U*P0*V; P1=U*P1*V; P2=U*P2*V;

C1=blkdiag(P2,I);
C0=[-P1/2 P2/2-P0/2;I O];
 



% % now I do one example as a sanity check
% QQ=orth(randn(3));
%  
% CC1=double(subs(C1,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{0 0 exp(-sqrt(42)) QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) 1}));
% CC0=double(subs(C0,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{0 0 exp(-sqrt(42)) QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) 1}));
%  
% MM=subs(M,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{0 0 exp(-sqrt(42)) QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) 1});
% vvf=subs(vf,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{0 0 exp(-sqrt(42)) QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) 1});
% s=solve([vvf(1) vvf(2) vvf(3)],{x y z});
% sz=sort(double(s.z));
% ez=sort(eig(CC0,CC1));
% norm(sz-ez)

 
CC1=vpa(subs(C1,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{x0_double y0_double z0_double QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) u_input}));
CC0=vpa(subs(C0,{x0 y0 z0 q11 q12 q13 q21 q22 q23 q31 q32 q33 u},{x0_double y0_double z0_double QQ(1,1) QQ(1,2) QQ(1,3) QQ(2,1) QQ(2,2) QQ(2,3) QQ(3,1) QQ(3,2) QQ(3,3) u_input}));

% CC0=U*CC0*V; CC1=U*CC1*V;

ez=(eig(CC0/CC1));

end
