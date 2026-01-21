function ind = get_ind_band(M,s)
    % This fonction return the linear index corresponding to the band of 
    % 2s+1 diag for a squared matrix (M by M). Thoses indexes are then 
    % stacked in a vector in the following order [d0; d1 ; d-1; d2 ....]
    %
    %----------------------------------------------------------
    % INPUT
    % M : size of squared matrix 
    % s : number of sub-diagonal and upper-diagonal to extract
    %----------------------------------------------------------
    % OUTPUT 
    % ind : vector containing linear index of K=2s+1 diagonales 
    %       stacked in the following order [d0; d1 ; d-1; d2, d-2 ....]
    
    [i, j] = meshgrid(1:M, 1:M);
    ind_mat = sub2ind([M, M], i', j');
    
    ind = diag(ind_mat,0);

    for k=(1:s)
        ind = [ind ; diag(ind_mat,k)];
        ind = [ind ; diag(ind_mat,-k)];
    end
    
end 