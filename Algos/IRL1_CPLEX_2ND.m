function[X, P_est, W, Fmin] = IRL1_CPLEX_2ND(Y, A, C, W0, X0, P0, par_problem, par_cplex, K, p)
    
    % Solve 2ND order problem with non-convexe lp-norm
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % C : OT cost matrix
    % W0, X0, P0 : initialization of variables & weights 
    % params_problem : hyperparameters of the problem 
    %	omega : OT term 
    % 	lambda : lp term
    % 	sdiag : bande parameter
    % par_solver : n° of algo used by CPLEX 
    %   barrier : 4
    % K : nbr of iterations for IRl1 scheme
    % p : value of p for lp-norm  
    %----------------------------------------------------------
    % OUTPUT
    %
    % X : solution (sequence of decompositions)
    % P_est : solution (set of OT plan)
    % W : final set of weights
    % Fmin : vectors of value of objective function 
    
    X{1} = X0;
    W{1} = W0; % weights initialization : 1st iteration eq to l1 norm
    
    Fmin = zeros(1,K);

    for k=1:K-1
        [X{k+1}, P_est, Fmin(k), ~, ~] = Algo_CPLEX_2ND_QP(Y, A, ...
            C, W{k}, X{k}, P0, par_problem, par_cplex); % Solve QP problem
        W{k+1} = p*X{k+1}.^(p-1); % Update weights
    end
end
