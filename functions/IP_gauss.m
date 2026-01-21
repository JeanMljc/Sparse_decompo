function h = IP_gauss(bin, mu, sigma, d)
    % Build Gaussian impulse respond (IP)
    %----------------------------------------------------------
    % INPUT 
    %
    % bin : grid
    % mu : mean value
    % sigma : standard deviation 
    % d (optional) : Gaussian is truncated beyond + or - d on the grid
    %----------------------------------------------------------
    % OUTPUT 
    %
    % h : gaussian IP, on grid bin
    
    h = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((bin-mu)/sigma).^2); % the IR is always centered 

    if nargin > 3
        h(abs(bin - mu) > d) = 0;
    end
end