function error = rMSE(x,x_star)
    % Compute root Mean Square Error
    %----------------------------------------------------------
    % INPUT 
    % x : solution
    % x_star : ground truth 
    %----------------------------------------------------------
    % OUTPUT 
    % error = norm(x-xc,'fro') / norm(x,'fro'); 
    
    error = norm(x(:)-x_star(:)) / norm(x_star(:)); 
end