function [M1,M2,A,B] = OT_sinkhorn_par(CK, X, K, param)
    
    % Solve problem type [Janati19] : solve N-1 OT problems 
    % with Sinkhorn algorithm
    %----------------------------------------------------------
    % INPUT 
    %
    % CK : repeated kernel corresponding to OT cost matrix
    % X : marginales of OT problem
    % X0 : intialization of variable X
    % param : hyperparameters of the problem 
    %   epsi : entropy term
    %   gamma : marginal relaxed term
    % K : numbers of iterations
    %----------------------------------------------------------
    % OUTPUT
    %
    % X2 : last iterate
    % dX : RMSE between consecutives iterates 
    
    [M,N] = size(X);
    
    X1 = X(:,1:N-1);
    X2 = X(:,2:N);

    B = ones(M,N-1);
    y = param.gamma / (param.epsi + param.gamma);
    
    for k=1:K
        A = (X1 ./ sum_MV(CK, B, 2)).^y;
        B = (X2 ./ sum_MV(CK, A, 1)).^y;
    end
    
    M1 = A .* sum_MV(CK, B, 2);
    M2 = B .* sum_MV(CK, A, 1);
end