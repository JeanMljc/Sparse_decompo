function X_re = recover_super(amp, supp, mu)
    
    % Build ground truth sequence of decompositions from dictionary grid 
    % and continuous trajectories
    %----------------------------------------------------------
    % INPUT
    % amp : amplitudes of the trajectories 
    % supp : continuous trajectories 
    % mu : dictionary grid 
    %----------------------------------------------------------
    % OUTPUT 
    % X_re : Ground truth decomposition on dictionary associated 
    % with mu
    
    [N,~] = size(amp);
    M = length(mu);
    
    % approximate continuous trajectories on discrete grid with
    % perfect_decompo()
    supp_grid = supp;
    supp_grid(:,2) = perfect_decompo(mu, supp(:,2));
    
    X_re = full(sparse(supp_grid(:,2), supp_grid(:,1), amp(:),M ,N));
end 
