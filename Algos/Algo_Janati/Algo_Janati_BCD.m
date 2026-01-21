function[X_est, P_est, dX] = Algo_Janati_BCD(Y, A, Ck, params_problem, param_algo)

    % Solve problem type [Janati19] : Entropy term and unbalanced term 
    % were added to criterion. Optimize with Block coordinate descente 
    % (BCD) algorithm 
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % Ck : kernel corresponding to OT cost matrix
    % params_problem : hyperparameters of the problem 
    %   omega : OT term 
    %   epsi : entropy term
    %   gamma : marginal relaxed term
    %   lambda : \ell_1 term 
    % par_solver : hyperparameters of algorithm
    %   K : number of iteration for the main algorithm
    %   K_prox : number of iteration for prox descente 
    %   K_sink : number of iteration OT algorithm
    %----------------------------------------------------------
    % OUTPUT
    %
    % X_est : solution (sequence of decompositions)
    % P_est : solution (set of OT plan)
    % dX : RMSE between consecutives iterates  

    [~,N] = size(Y);
    [~,M] = size(A);

    CK = repmat(Ck, 1, 1, N-1);

    % Get hyperparameters of algorithm
    K = param_algo.K;
    K_prox = param_algo.K_prox;
    K_sink = param_algo.K_sink;
    
    % Initialization

    X = cell(K,1);
    dX = cell(K,1);

    X{1} = ones(M,N);
    M1 = ones(M,N-1);
    M2 = ones(M,N-1);
    
    for k = 1:K
        [X{k+1}, ~] = descent_prox(Y, A, X{k}, M1, M2, params_problem, K_prox);
        [M1,M2,u,v] = OT_sinkhorn_par(CK, X{k+1}, K_sink, params_problem);
        dX{k} = norm(X{k+1}-X{k}, 'fro') / norm(X{k+1}, 'fro');
    end
    
    % Output variables

    X_est = X{end};
    dX = cell2mat(dX);

    P_est = zeros(M,M,N-1);
    for n=1:N-1
        P_est(:,:,n) = repmat(u(:,n),1,M).*Ck.*repmat(v(:,n)',M,1);
    end
    
end 