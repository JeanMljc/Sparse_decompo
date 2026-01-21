function [X_sol,V_sol,Phi] = Algo_dynamic_Elvander(Y, A, CostMat, NV, param)
    
    % Script from paper [Elvander20] 
    % This corresponds to Example 1 and Section 6.1 of
    %
    % "Multi-marginal optimal transport using partial information with
    % applications in robust localization and sensor fusion", 2020, Elvander et
    % al.
    
    NTheta = size(A,2);
    N = size(Y,2);
        
    %% Dimension reduction of data and operators
    
    r_cell_real = cell(1,N);
    for n=1:N
        r_cell_real{n} = Y(:,n);
    end
    
    %% Compute covariance operators
    
    % if we introduce a hidden velocity state, the covariance operators must
    % include marginalization of the velocity spectra. Also the adjoint and
    % Jacobian must reflect this. More efficient alternative to matrix
    % representation.
    vec_op = @(var) var(:);
    operator_cell = cell(1,length(r_cell_real));
    adjoint_cell = cell(1,length(r_cell_real));
    Jacobian_cell = cell(1,length(r_cell_real));
    
    for kcell = 1:length(r_cell_real)
        G_operator = @(var) A*sum(reshape(var(:),NV,NTheta))';
        operator_cell{kcell} = G_operator;
        G_adjoint = @(var) reshape(repmat(vec_op(A'*var)',NV,1),NV*NTheta,1);
        adjoint_cell{kcell} = G_adjoint;
        G_Jacobian = @(var) A*diag(sum(reshape(var(:),NV,NTheta)))*A';
        Jacobian_cell{kcell} = G_Jacobian;
    end
        
    %% Solve OMT tracking problem with dynamics
    
    % Entropy regularization parameter
    epsilon = param.epsilon;
  
    % Regularization parameter penalizing measurement error
    gamma = param.gamma;
    
    % Solver tolerance
    tol = param.tol;
    
    % CostMat = CostMat ./ 100; % ATTENTION AJOUT JEAN M
    
    Cost_cell = cell(1); Cost_cell{1} = CostMat;
    fprintf('\n')
    fprintf('----------------------------------------------\n')
    fprintf('Solving OMT tracking problem with dynamics...\n')
    [Phi,~,~] = ...
        multimarginal_tracking_sinkhornnewton(operator_cell,adjoint_cell,Jacobian_cell,r_cell_real,Cost_cell, epsilon,gamma,tol);
    fprintf('----------------------------------------------\n')
    
    %% Plot estimates
    
    N_Phi = size(Phi,1);
    Kv = N_Phi/NV;
    
    X_sol = zeros(NTheta,N);
    V_sol = zeros(NV,N);
    
    for n=1:N
        phi_n = Phi(:,n);
        x_sol = zeros(NTheta,1);
        v = zeros(NV,1);
        for k=1:Kv
            ph = phi_n((k-1)*NV+1:k*NV);
            v = v + ph;
            x_sol(k) = sum(ph);
        end
        V_sol(:,n) = v;
        X_sol(:,n) = x_sol;
    end
end