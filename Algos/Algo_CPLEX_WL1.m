function [X_est, F_cplex, exitflag, output] = Algo_CPLEX_WL1(Y, A, W, X0, params, par_solver)

    % Solve BPDN problem : QP problem
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % W : weights for L1 norm (with W = ones(M,N) => BPDN/LASSO)
    % X0 : variables initialization
    % params_problem : hyperparameters of the problem 
    % 	lambda : l1 term
    % par_solver : n° of algo used by CPLEX 
    %   auto : 0
    %   primal : 1
    %   primal : 2
    %   barrier : 4
    %----------------------------------------------------------
    % OUTPUT
    %
    % X_est : solution (sequence of decompositions)
    % F_cplex : minimum of objective value
    % exitflag : exitflag of CPLEX solveur
    % output : output of CPLEX solveur

    lambda = params.lambda; 

    [~,N] = size(Y);
    [~,M] = size(A);
    z = N*M;
    
    % Initialization
    z_init = X0(:);

    % QP standard CPLEX form : min  0.5*x'*Q*x + q*x

    % Q : quadratic term

    Q = 2*kron(speye(N),A'*A);

    % q : linear term

    q = lambda*W - 2*(A')*Y;
    q = q(:); 

    % Positivity constraint : lb (lower bound)
    
    lb = zeros(z,1);

    % Calling CPLEX %
        
    options = cplexoptimset('cplex');
    options.display = par_solver.display;
    options.qpmethod = par_solver.algo;
    
    [v,fval,exitflag,output]= cplexqp(Q, q', [], [], [], [], lb, [], z_init, options);

    % Solution : decompositions X_est
    X_est = reshape(v,M,N);

    % Value of objective function
    F_cplex = fval + Y(:)'*Y(:); % add constant term Y'Y

end
