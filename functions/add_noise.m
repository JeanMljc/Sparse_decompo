function [Yn, B, sigma_noise] = add_noise(Y,SNR)
    % Build noisy signal Yn from ground truth Y for given SNR
    %----------------------------------------------------------
    % INPUT 
    % Y : ground truth signal 
    % SNR : required signal-to-noise ratio (in dB)
    %----------------------------------------------------------
    % OUTPUT 
    %
    % Yn : noisy signal with Gaussian noise 
    % B : Gaussian noise term 
    % sigma_noise : standard deviation of the noise
    
    N = numel(Y);
    Py = (1/N)*sum(Y(:).^2);
    sigma_noise = sqrt(Py*10^(-SNR/10));
    B = sigma_noise*randn(size(Y));
    Yn = Y + B;
end