function X_re = recover(amp, supp, M, N, d)
    % Build groud truth decompositions on the dictionary using continuous
    % trajectories & amplitude 
    %
    % WARNING : DO NOT WORK WITH OVERSAMPLED GRID
    %----------------------------------------------------------
    % INPUT 
    %
    % amp : ground amplitude
    % supp : ground truth support (continuous trajectories)
    % M, N : size & number of decompositions 
    % d (optional) : Gaussian is truncated beyond + or - d on the grid
    %----------------------------------------------------------
    % OUTPUT 
    %
    % X_re : Ground truth sequence of N decompositions (M,N)
    
    % approximate trajectory on dictionary grid
    supp_grid = supp;
    supp_grid(:,2) = round(supp(:,2))-d;
    
    % build corresponding decompositions
   
    X_re = full(sparse(supp_grid(:,2), supp_grid(:,1), amp(:) ,M, N));
end 
