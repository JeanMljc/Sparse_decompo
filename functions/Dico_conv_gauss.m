function [A, mu, d_tail] = Dico_conv_gauss(N_lambda, sigma, kappa)

    % Build a convolution gaussian dictionary from a Gaussian IP determine
    % by sigma.
    %
    % Hypothesis "full" for the convolution
    %----------------------------------------------------------
    % INPUT 
    % 
    % N_lambda : numbers of points on observation grid 
    % sigma : standard deviation of of Gaussian IP
    % kappa : over-sampling factor for the dictionary grid
    % kappa : Gaussian is truncated beyond + or - d on the grid
    %----------------------------------------------------------
    % OUTPUT 
    %
    % A : convolution gaussian dictionary of size (N_lambda,M)
    % mu : grid of dictionary A (M,1)
    % d_tail : the gaussian is truncated beyond ± d_tail
    
    % define step of the dictionary grid
    Delta_mu = 1 / kappa;

    % build dictionary grid
    d_tail = floor(3*sigma);
    mu = d_tail+1:Delta_mu:N_lambda-d_tail;

    % build dictionary
    M = length(mu);
    A = zeros(N_lambda,M);

    for m=1:M
        A(:,m) = IP_gauss(1:N_lambda, mu(m), sigma, d_tail);
    end
    
end 