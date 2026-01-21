function[x_est, P_est, F_cplex, exitflag, output] = Algo_CPLEX_1ST_QP(Y, A, C, W, X0, P0, params_problem, par_solver)
    
    % Solve 1ST order problem : QP problem
    % NO MORE SUM-UP TO ONE ASSUMPTION 
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % C : OT cost matrix
    % W : weight vector for l1 norm
    % X0, P0 : initialization of variables 
    % params_problem : hyperparameters of the problem 
    %	omega : OT term 
    % 	lambda : l1 term
    % 	sdiag : bande parameter
    % par_solver : n° of algo used by CPLEX 
    %   auto : 0
    %   primal : 1
    %   primal : 2
    %   barrier : 4
    %----------------------------------------------------------
    % OUTPUT
    %
    % x_est : solution (sequence of decompositions)
    % P_est : solution (set of OT plan)
    % F_cplex : minimum of objective value
    % exitflag : exitflag of CPLEX solveur
    % output : output of CPLEX solveur
    
    omega = params_problem.omega;
    lambda = params_problem.lambda;
    sdiag = params_problem.sdiag;
    
    [~,N] = size(Y);
    [~,M] = size(A);

    z_init = [X0(:) ; P0(:)]; % WARNING DO NOT SUPPORT BAND STRUCTURE sdiag<M-1
    
    % QP standard CPLEX form : min  0.5*z'*Q*z + q*z sc Aeq*z = beq, z >= 0

    zx = N*M; % size of x
    zp = M + 2*sdiag*M - sdiag*(sdiag+1);  % size of P_band
    z = zx + (N-1)*zp; % size of variable z

    % Get indexes of the band
    ind_band = get_ind_band(M,sdiag);

    % Band cost Matrix for OT term 
    
    cost = C(ind_band);
    
    % Q : quadratic term

    Q = sparse(z,z);
    Q(1:zx,1:zx) = 2*kron(speye(N),A'*A);

    % q : linear term

    LD = lambda*W - 2*A'*Y;
    LC = sparse(repmat(omega.*cost,N-1,1));
    q = [LD(:);LC];

    clear("LD","LC")

    % Equality constraints : Aeq & beq 
        
    M1 = [sparse((N-1)*M,M), speye((N-1)*M)];
    M2 = [speye((N-1)*M), sparse((N-1)*M,M)];
    
    s1 = kron(speye(M),ones(1,M));
    s2 = kron(ones(1,M),speye(M));

    % Adapt operator to band storage

    S1 = kron(speye(N-1),s1(:,ind_band));
    S2 = kron(speye(N-1),s2(:,ind_band));
    
    Aeq = [ -M2, S2 ; -M1, S1 ];
    beq = sparse(2*M*(N-1),1);
    
    clear("M1","M2","S1","S2","s1","s2")

    % Positivity constraint : lb (lower bound)

    lb = sparse(z,1);

    % Calling CPLEX %

    options = cplexoptimset('cplex');
    options.display = par_solver.display;
    options.qpmethod = par_solver.algo;
        
    [v, fval, exitflag, output]= cplexqp(Q, q', [], [], Aeq, beq ,lb, [], z_init, options);
        
    % Solution : decompositions X_est
    x_est = reshape(v(1:zx),M,N);
    
    % Solution : OT plans P_est
    vp = v(zx+1:end);
    P_est = zeros(M,M,N-1);
    
    for n=1:N-1
        Pn = zeros(M,M);
        Pn(ind_band) = vp((n-1)*zp+1:n*zp);
        P_est(:,:,n) = Pn;
    end

    % Value of objective function    
    F_cplex = fval + Y(:)'*Y(:); % add constant term Y'Y
    
end
