function X_sol = Algo_static_Elvander(Y, A, CostMat, param)
    
    % Script from paper [Elvander20] 
    % This corresponds to Example 1 and Section 6.1 of
    %
    % "Multi-marginal optimal transport using partial information with
    % applications in robust localization and sensor fusion", 2020, Elvander et
    % al.
    
    % NTheta = size(A,2);
    N = size(Y,2);
        
    %% Dimension reduction of data and operators
    
    r_cell_real = cell(1,N);
    for n=1:N
        r_cell_real{n} = Y(:,n);
    end
    
    %% Solve OMT problem without dynamics (no velocity state)

    % We only have the postion state: covariance operator and related operators
    % can be directly represented as matrices.
    operator_cell_static = cell(1,length(r_cell_real));
    adjoint_cell_static = cell(1,length(r_cell_real));
    Jacobian_cell_static = cell(1,length(r_cell_real));
    for kcell = 1:length(r_cell_real)
        G_operator = @(mu) A*mu;
        operator_cell_static{kcell} = G_operator;
        G_adjoint = @(mu) A'*mu;
        adjoint_cell_static{kcell} = G_adjoint;
        G_Jacobian = @(mu) A*diag(mu)*A';
        Jacobian_cell_static{kcell} = G_Jacobian;
    end
    
    Cost_cell_static = cell(1,1);
    Cost_cell_static{1} = CostMat;

    % Entropy regularization parameter
    epsilon = param.epsilon;
  
    % Regularization parameter penalizing measurement error
    gamma = param.gamma;
    
    % Solver tolerance
    tol = param.tol;
        
    fprintf('\n')
    fprintf('----------------------------------------------\n')
    fprintf('Solving OMT tracking problem without dynamics...\n')
    [X_sol,~,~] = ...
        multimarginal_tracking_sinkhornnewton(operator_cell_static,adjoint_cell_static,...
    Jacobian_cell_static,r_cell_real,Cost_cell_static, epsilon, gamma, tol);
    fprintf('----------------------------------------------\n')
    
%     %% Plot estimates
%     
%     % Plot estimate using dynamic model as well as static model
%     % nbrSnapshots = N;
%     % plot_estimate
%     
%     N_Phi = size(Phi,1);
%     Kv = N_Phi/NV;
%     
%     X_sol = zeros(NTheta,N);
%     V_sol = zeros(NV,N);
%     
%     for n=1:N
%         phi_n = Phi(:,n);
%         x_sol = zeros(NTheta,1);
%         v = zeros(NV,1);
%         for k=1:Kv
%             ph = phi_n((k-1)*NV+1:k*NV);
%             v = v + ph;
%             x_sol(k) = sum(ph);
%         end
%         V_sol(:,n) = v;
%         X_sol(:,n) = x_sol;
%     end
end