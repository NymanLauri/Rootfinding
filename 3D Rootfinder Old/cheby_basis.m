function y = cheby_basis(n)
y= zeros(n);
y(n,1) = 1;
y(n-1,2) = 1;

for i = 3:n
   y(1:end-1,i) = 2*y(2:end,i-1);
   y(:,i) =  y(:,i) - y(:,i-2);
end

end