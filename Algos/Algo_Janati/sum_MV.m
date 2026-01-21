function res = sum_MV(T, A, mode)

    % Perform matrix-vector product : each slice T(:,:,n) 
    % with each column vector A(:,n).
    %----------------------------------------------------------
    % INPUT 
    %
    % T : Tensor containing matrix/slices (M, K, N)
    % A : Matrice containing vector (K, N) or (M, N)
    % mode : mode of slice to sum on (1 : T(:,:,n)^{\top} or 2 : T(:,:,n)) 
    %----------------------------------------------------------
    % OUTPUT
    %
    % res : Matrix obtained (K, N) or (M, N)
    
    [M,N] = size(A);
    res = zeros(M,N);
    
    % Get the right mode to multiply
    if mode == 1
        T = permute(T,[2 1 3]);
    end 
    % matrix-vector product
    for n=1:N
        res(:,n) = T(:,:,n)*A(:,n);
    end
end