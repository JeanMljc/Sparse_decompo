function[T, P, X, F_cplex, exitflag, options] = Algo_interpo_R2(C, Marg, X0, N, par_cplex)
    
    % Solve interpolation problem with R2 : LP problem.
    % Working on full OT variables (no thin storage here) 
    %----------------------------------------------------------
    % INPUT 
    %
    % C : OT cost tensor (M x M x M)
    % Marg : first & last decomposition (M x 2)
    % X0 : initialization of variable X
    % N : number of time step over which interpolation is made
    % par_solver : n° of algo used by CPLEX 
    %   auto : 0
    %   primal : 1
    %   dual : 2
    %   barrier : 4
    %----------------------------------------------------------
    % OUTPUT
    %
    % T : sequence of OT tensors
    % P : sequence of OT matrix (OT plan)
    % X : sequence of decompositions
    % F_cplex : minimum of objective value
    % exitflag : exitflag of CPLEX solveur
    % output : output of CPLEX solveur
        
    M = size(Marg,1);
    
    z_x = M*N;
    z_P = (N-1)*M^2;
    z_T = (N-2)*M^3;
    z = z_x + z_P + z_T;

    % Initialization
    z_init = [X0(:) ; zeros(z_P + z_T,1)]; 

    C_rep = repmat(C,1,1,1,N-2);

    % Objective function

    q = [sparse(z_x,1) ; sparse(z_P,1) ; C_rep(:)];

    % Equality constraints : interpolation
    
    % Linear operator L1 select decompostion x_1 and x_N

    L1a = sparse(2*M,z_x);
    L1a(1:M,1:M) = speye(M);
    L1a(M+1:end,end-M+1:end) = speye(M);
    L1 = sparse(2*M,z);
    L1(:,1:z_x) = L1a;

    clear("L1a")

    % Equality constraints OT : Build L0
    
    M1 = [sparse((N-1)*M,M), speye((N-1)*M)];
    M2 = [speye((N-1)*M), sparse((N-1)*M,M)];
    
    s1 = kron(speye(M),ones(1,M));
    s2 = kron(ones(1,M),speye(M));
    
    S1 = kron(speye(N-1),s1);
    S2 = kron(speye(N-1),s2);
    
    L01 = sparse(2*(N-1)*M,z);
    L01(:,1:z_x+z_P) = [ -M2, S2 ; -M1, S1 ];
    
    clear("M1","M2","S1","S2","s1","s2")

    % Equality constraints OT : Build L1
    
    K1 = [sparse((N-2)*M^2,M^2), speye((N-2)*M^2)];
    K3 = [speye((N-2)*M^2), sparse((N-2)*M^2,M^2)];

    sig1 = kron(speye(M^2),ones(1,M));
    sig3 = kron(ones(1,M),speye(M^2));
    
    Sig1 = kron(speye(N-2),sig1);
    Sig3 = kron(speye(N-2),sig3);

    L02 = sparse(2*(N-2)*M^2,z);
    L02(:,z_x+1:end) = [ -K3, Sig3 ; -K1, Sig1 ];

    clear("K1","K3","sig1","sig3","Sig1","Sig3")
    
    % Build Aeq and beq

    Aeq = [L1 ; L01 ; L02];

    beq = sparse(size(Aeq,1),1);
    beq(1:2*M) = Marg(:); % fix equality constraint on x_1 and x_N

    % Lower bound : positivity constraints

    lb = sparse(z,1);

    % Calling CPLEX %

    options = cplexoptimset;
    options.display = par_cplex.display;
    options.lpmethod = par_cplex.algo;

    [z_sol, F_cplex, exitflag, options] = cplexlp(q, [], [], Aeq, beq ,lb, [], z_init, options);

    % test cplexlp  instead of cplexqp
    
    X_est = z_sol(1:z_x);
    P_est = z_sol(z_x+1:z_x+z_P);
    T_est = z_sol(z_x+z_P+1:end);

    % Solution X_cplex
    X = reshape(X_est,M,N);
    P = reshape(P_est,M,M,N-1);
    T = reshape(T_est,M,M,M,N-2);
end