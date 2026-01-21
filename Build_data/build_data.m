%% This script need theses variables to run properly :

% SNR 
% N 
% N_lambda 
% std_problem 

%% Build continuous trajectories

% time axis
times_L = linspace(0,1,N);
times_Q = linspace(-1,1,N);

Poly = cell(1,1);

% Linear trajectories
Poly{1} = [7,24]; 
Poly{2} = [-10,32];

% Quadratic trajectories
Poly{3} = [6,5,18];

% Get number of trajectories
N_tra = length(Poly);

%% Build discret trajectories

tra = zeros(N,N_tra);

for i=1:N_tra
    coeff = Poly{i};
    if length(coeff)==2
        tra(:,i) = coeff(1)*times_L + coeff(2);
    else
        tra(:,i) = coeff(1)*times_Q.^2 + coeff(2)*times_Q + coeff(3);
    end 
end

% Define weight of each trajectory and total mass 
weights = [2 3 2];
mass = sum(weights);

% check consistency
if (length(Poly)~=N_tra)||(length(weights)~=N_tra)
    error("Not the same number of trajectories and weights/tra");
end

%% Build multichannel signal Z

% Z is the convolution of trajectories with Gaussian IP

% Observation grid
z_bin = 1:N_lambda;
K = length(z_bin);

% Build multichannel signal Z  
Z_star = zeros(K,N);

for n=1:N
    z_channel = zeros(K,1);
    for i=1:N_tra
        z_channel = z_channel + weights(i)*IP_gauss(z_bin, tra(n,i), std_problem)';
    end    
    Z_star(:,n) = z_channel;
end

%% Build groud truth continuous trajectories

supp_star = [repmat((1:N)',N_tra,1),tra(:)];
amp_star = repmat(weights,N,1);

%% Plot

% figure()
% imagesc(Z_star')
% hold on
% plot(tra',1:N,'r')
% title("Data Z^*")
% axis equal
% axis tight
% colorbar()
