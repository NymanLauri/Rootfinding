% Computes the division between the polynomial B in variables s1,s2,t1,t2,z and the term (s1-t1)(s2-t2).
% Polynomial B is given in tensor form. Returns the resulting tensor.
function quotient = cayley_divide(B)

    [n1, n2, n3, n4, n5] = size(B);
    assert(n1==n3)
    assert(n2==n4)
    
    % We know that B is divisible by x-z, so then the degree of z cannot at any point exceed
    % nx when performing long division (dividing in variable x). Thus, no
    % need for larger tensors and (nx,ny,nz) suffices.
    quotient = zeros(n1,n2,n3,n4,n5);
    temp = B;
    
    % Perform the division
    for i = flip(2:n1)
        quotient(i-1,:,:,:,:) = temp(i,:,:,:,:);
        temp(i-1,:,2:end,:,:) = temp(i-1,:,2:end,:,:)+quotient(i-1,:,1:end-1,:,:);
        temp(i,:,:,:,:) = temp(i,:,:,:,:)-quotient(i-1,:,:,:,:); % This line could be left out
    end
 
    temp = quotient;
    for i = flip(2:n2)
        quotient(:,i-1,:,:,:) = temp(:,i,:,:,:);
        temp(:,i-1,:,2:end,:) = temp(:,i-1,:,2:end,:)+quotient(:,i-1,:,1:end-1,:);
        temp(:,i,:,:,:) = temp(:,i,:,:,:)-quotient(:,i-1,:,:,:); % This line could be left out
    end
    
    % Maximal degrees of s1,s2,t1 and t2 get decreased by 1 during the division
    % so the resulting zero coefficients can be left out
    quotient = quotient(1:end-1,1:end-1,1:end-1,1:end-1,:);
end

