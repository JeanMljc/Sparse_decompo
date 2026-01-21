function[X2, dX] = descent_prox(Y, A, X0, M1, M2, param, kmax)

    % Solve problem type [Janati19] : perform proximal descent
    % with algo [Fercoq15]
    %----------------------------------------------------------
    % INPUT 
    %
    % Y : signal to decompose
    % A : sparse dictionary
    % X0 : intialization of variable X
    % param : hyperparameters of the problem 
    %   omega : OT term 
    %   epsi : entropy term
    %   gamma : marginal relaxed term
    %   lambda : \ell_1 term 
    % k_max : numbers of iterations for proximal descente
    %----------------------------------------------------------
    % OUTPUT
    %
    % X2 : last iterate
    % dX : RMSE between consecutives iterates  
    
    N = size(Y,2);
    M = size(A,2);

    % Get hyperparameters of problem 

    par_g = param.omega*param.gamma;
    lamb = param.lambda;

    b = 2*par_g*ones(M,1) + lamb;
    b(1) = par_g + lamb;
    b(M) = par_g + lamb;
    
    % Get hyperparameters of algorithm
    tau = 1;
    % tau = 1 / sum(A,1);
    
    % Initialization of algorithm

    Mar = par_g*([M1,zeros(M,1)] + [zeros(M,1),M2]);
    dX = cell(10,1);
    X1 = X0;
    k = 1;
    
    % Proximal descente
    
    while true
        
        X2 = zeros(M,N);

        for n=1:N
            X1v = X1(:,n);
            y = Y(:,n);
            a = Mar(:,n);
            for m=1:M
                x0 = X1v(m) - tau*dot(A(:,m),(A*X1v - y));
                x2 = 0.5*(x0-b(n) + ((x0-b(n))^2 + 4*a(m))^(0.5));
                X2(m,n) = x2;
            end
        end
        
        dX{k} = norm(X2-X1, 'fro') / norm(X2, 'fro'); 
        % F{k} = norm(Y-A*X2,"fro")^2 + eval_f(Mar,X2,par_g);

        if k > kmax
            dX = cell2mat(dX);
            break
        end
        
        X1 = X2;
        k = k + 1; 
    end

end 