function[X, W, Fmin] = IRL1_CPLEX_LP(Y, A, W0, X0, par_problem, par_cplex, K, p)

    % Solve lp regularized problem with non-convexe lp-norm
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % W0, X0 : initialization of variables & weights 
    % params_problem : hyperparameters of the problem 
    % 	lambda : lp term
    % par_solver : n° of algo used by CPLEX 
    %   barrier : 4
    % K : nbr of iterations for IRl1 scheme
    % p : value of p for lp-norm  
    %----------------------------------------------------------
    % OUTPUT
    %
    % X : solution (sequence of decompositions)
    % W : final set of weights
    % Fmin : vectors of value of objective function  
        
    X{1} = X0;
    W{1} = W0; % weights initialization : 1st iteration eq to l1 norm
    
    Fmin = zeros(1,K);

    for k=1:K-1
        [X{k+1}, Fmin(k), ~, ~] = Algo_CPLEX_WL1(Y, A, W{k} ...
            , X{k}, par_problem, par_cplex); % Solve QP problem
        W{k+1} = p*X{k+1}.^(p-1); % Update weights
    end
end
