function [x_est, P_est, F_cplex, exitflag, output] = Algo_CPLEX_2ND_QP(Y, A, C, W, X0, P0, params, par_solver)
    
    % Solve 2ND order problem : QP problem
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
    
    omega = params.omega;
    lambda = params.lambda;
    idx_3D = params.ind;

    [~,N] = size(Y); 
    [~,M] = size(A);

    % Initialization
    z_init = [X0(:) ; P0(:)]; 
    
    % QP standard CPLEX form : min  0.5*z'*Q*z + q*z sc Aeq*z = beq, z > 0
    
    zx = N*M; % size of x
    zp = (N-1)*M^2; % size of 2D coupling : full
    zP = (N-2)*numel(C); % size of 3D coupling : band
    
    z = zx + zp + zP;
    
    % WARNING : OT plans (2D coupling) is not store in band !
    
    % Q : quadratic term

    Q = sparse(z,z);
    Q(1:zx,1:zx) = 2*kron(speye(N),A'*A);

    % L : linear term

    LD = lambda*W - 2*A'*Y;
    LC = sparse(repmat(omega.*C,N-2,1));
    L = [LD(:); sparse(zp,1); LC];

    clear("LC","LD")
    
    % Equality constraints : Build L0
    
    M1 = [sparse((N-1)*M,M), speye((N-1)*M)];
    M2 = [speye((N-1)*M), sparse((N-1)*M,M)];
    
    s1 = kron(speye(M),ones(1,M));
    s2 = kron(ones(1,M),speye(M));
    
    S1 = kron(speye(N-1),s1);
    S2 = kron(speye(N-1),s2);
    
    L1 = sparse(2*(N-1)*M,z);
    L1(:,1:zx+zp) = [ -M2, S2 ; -M1, S1 ];
    
    clear("M1","M2","S1","S2","s1","s2")
    
    % Equality constraints : Build L1
    
    K1 = [sparse((N-2)*M^2,M^2), speye((N-2)*M^2)];
    K3 = [speye((N-2)*M^2), sparse((N-2)*M^2,M^2)];

    sig1 = kron(speye(M^2),ones(1,M));
    sig3 = kron(ones(1,M),speye(M^2));
    
    Sig1 = kron(speye(N-2),sig1(:,idx_3D));
    Sig3 = kron(speye(N-2),sig3(:,idx_3D));
    
    L2 = sparse(2*(N-2)*M^2,z);
    L2(:,zx+1:end) = [ -K3, Sig3 ; -K1, Sig1 ];

    clear("K1","K3","sig1","sig3","Sig1","Sig3")
    
    % Equality constraints : Build Aeq and beq

    Aeq = [L1 ; L2];
    beq = sparse(size(Aeq,1),1);
    
    % Positivity constraint : lb (lower bound)

    lb = sparse(z,1);  
    
    % Print size of problem
    
    fprintf('Size of global variable z = %d \n', z);
    
    % Calling CPLEX %

    options = cplexoptimset('cplex');
    options.display = par_solver.display;
    options.qpmethod = par_solver.algo;

    [v, fval, exitflag, output]= cplexqp(Q, L', [], [], Aeq, beq ,lb, [], z_init, options);

    % Solution : decompositions x_est
    x_est = reshape(v(1:zx),M,N);
    
    % Solution : OT plans P_est
    vp = v(zx+1:end);
    
    P_est = zeros(M,M,N-1);
    for n=1:N-1
        P_est(:,:,n) = reshape(vp((n-1)*M^2+1:n*M^2),M,M);
    end
    
    % Value of objective function
    F_cplex = fval + Y(:)'*Y(:); % add constant term Y'Y

end
