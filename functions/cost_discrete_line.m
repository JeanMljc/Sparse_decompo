function [f, cost] = cost_discrete_line(D,N)

    % Compute R2 cost of discrete line numerically 
    %----------------------------------------------------------
    % INPUT
    % D : drift between first and last decom
    % N : Number of decompositions to interpolate
    %----------------------------------------------------------
    % OUTPUT 
    % f : vector of cost at each channel n 
    % cost : total average cost 
    
    alpha = D/N;
    f = zeros(N,1);
    
    for n=1:N
        f(n) = abs(2*floor(alpha*(n+1)) -floor(alpha*n) -floor(alpha*(n+2)));
    end
    
    cost = (1/N)*sum(f);
end
