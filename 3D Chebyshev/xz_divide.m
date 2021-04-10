% Computes the division between the polynomial B in variables x,y,z and the term x-z.
% Polynomial B is given in tensor form. Returns the resulting tensor.
function quotient = xz_divide(B)

    [nx, ny, nz] = size(B);
    assert(nx==nz)
    
    % We know that B is divisible by x-z, so then the degree of z cannot at any point exceed
    % nx when performing long division (dividing in variable x). Thus, no
    % need for larger tensors and (nx,ny,nz) suffices.
    quotient = zeros(nx,ny,nz);
    temp = B;
    
    % Perform the division
    for i = flip(2:nx)
        quotient(i-1,:,:) = temp(i,:,:);
        temp(i-1,:,2:end) = temp(i-1,:,2:end)+quotient(i-1,:,1:end-1);
        temp(i,:,:) = temp(i,:,:)-quotient(i-1,:,:); % This line could be left out
    end
    
    % Maximal degrees of x and z get decreased by 1 during the division
    % so the resulting zero coefficients can be left out
    quotient = quotient(1:end-1,:,1:end-1);
end

