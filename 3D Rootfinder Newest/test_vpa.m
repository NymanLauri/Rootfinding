% Error behaviour of System 1 from "bad systems" with different arithmetic
% precisions
interval = [10.^(-32:0)];
case1 = zeros(1,length(interval));
case2 = zeros(1,length(interval));
case3 = zeros(1,length(interval));
for k=1:50
k
vals1 = [];
vals2 = [];
vals3 = [];
% QQ=orth(randn(3));
U=orth(randn(4)); V=orth(randn(4));
% U=eye(4); V=eye(4);
for u = interval
  
    x0 = 2*rand-1;
    y0 = 2*rand-1;
    z0 = 2*rand-1;

    % Orthogonal matrix 1
    QQ=[sqrt(1/3) sqrt(1/3) sqrt(1/3); sqrt(2/3) -sqrt(1/6) -sqrt(1/6); 0 sqrt(1/2) -sqrt(1/2)];
%     QQ=orth(randn(3));
    
    myflag = true;
    while myflag
        try
            digits(16)
            roots1 = exact_resultant(QQ, x0, y0, z0, u, U, V);

            digits(32)
            roots2 = exact_resultant(QQ, x0, y0, z0, u, U, V);

            digits(64)
            roots3 = exact_resultant(QQ, x0, y0, z0, u, U, V);    
            
            vals1 = [vals1 min(abs(roots1 - z0))];
            vals2 = [vals2 min(abs(roots2 - z0))];
            vals3 = [vals3 min(abs(roots3 - z0))];
            myflag=false;
        catch
            fprintf('No convergence in iteration k = %d with u = %d, redoing.\n', k, u);
            x0 = 2*rand-1;
            y0 = 2*rand-1;
            z0 = 2*rand-1;
        end
    end
end

case1 = case1 + vals1/k;
case2 = case2 + vals2/k;
case3 = case3 + vals3/k;

end

figure;
loglog(interval,case1, '-x',interval,case2, '-x',interval,case3, '-x')
legend('digits=16','digits=32','digits=64')
xlabel('$u$','fontsize',16,'Interpreter','latex')
ylabel('Inaccuracy of the $z$-component','fontsize',16,'Interpreter','latex')
xlim([interval(1),interval(end)])
