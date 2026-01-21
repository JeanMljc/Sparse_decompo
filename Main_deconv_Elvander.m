clear
close all

config

%% Build ground truth multichannel data : Z*

% Problem setting : size of data, number of sources, Gaussian IP
N = 20;
N_lambda = 60; 
std_problem = 3; 

% Build observations Y
run("Build_data/build_data.m");   

%% Add noise to ground truth multichannel data : Z

% noise level
SNR = 30;

% number of noise realisation 
Nb_noise = 10;
Zn = cell(8,1);

rng(0);
for i=1:Nb_noise
    Zn{i} = add_noise(Z_star,SNR);
end

Z = Zn{5};

%% Build Dictionary : A

% define step of the dictionary grid
kappa = 1;

% Build convolutive Gaussian dictionary
[A, mu, d_tail] = Dico_conv_gauss(N_lambda, std_problem, kappa);
M = length(mu);

%% Build ground truth signal & decomposition : Z_GT, X_GT 

% Ground-truth decompositions on the grid (does not handle oversampled grids)
X_GT = recover(amp_star, supp_star, M, N, d_tail);
Z_GT = A*X_GT;

% Error between Z_GT and Z_star
rMSEz_GT = rMSE(Z_GT,Z_star);

%% Static method

% Grid of position
gridTheta = (1:M)'/ M; % why divide by M ?

% Cost matrix : static
ett_theta = ones(M,1);
Cost_static = abs(gridTheta*ett_theta'-ett_theta*gridTheta').^2;

% Hyperparameters static method

par_FE_static.epsilon = 10^-3;
par_FE_static.gamma = 1; % ATTENTION : 1/omega
par_FE_static.tol = 1e-3;

% Call static method [Elvander20]

X1_FE = Algo_static_Elvander(Z, A, Cost_static, par_FE_static);

Z1_FE = A*X1_FE;
rMSEx_FE1 = rMSE(X1_FE,X_GT);
rMSEz_FE1 = rMSE(Z1_FE,Z_star);    

% figure()
% imagesc(X1_FE')

%% Dynamic method

% Grid of speed
NV = 50; 
vRange = [-0.05,0.05];
gridV = linspace(vRange(1),vRange(2),NV);

% Cost matrix : dynamic

Winv = [12, -6; -6, 4];
eA1 = [1,1;0,1];

Cost_dynamic = zeros(M*NV,M*NV);

fprintf('Computing cost matrix for dynamical model...')
for kTheta1 = 1:M
    for kV1 = 1:NV
        index1 = (kTheta1-1)*NV + kV1;
        x1 = [gridTheta(kTheta1);gridV(kV1)];
        for kTheta2 = 1:M
           for kV2 = 1:NV
               index2 = (kTheta2-1)*NV + kV2;
               x2 = [gridTheta(kTheta2);gridV(kV2)];
               costVal = (x2-eA1*x1)'*Winv*(x2-eA1*x1);
               Cost_dynamic(index1,index2)= costVal;
           end
        end 
    end
end
fprintf('done.\n')

%% Hyperparameters dynamic method

par_FE.epsilon = 10^-4;
par_FE.gamma = 10^-2; % ATTENTION : 1/omega
par_FE.tol = 1e-3;

%% Call Elvander solver

[X2_FE, v, phi] = Algo_dynamic_Elvander(Z, A, Cost_dynamic, NV, par_FE);

Z2_FE = A*X2_FE;
rMSEx_FE2 = rMSE(X2_FE,X_GT);
rMSEz_FE2 = rMSE(Z2_FE,Z_star);

%% Plots section

figure(1)
imagesc(X1_FE')
title("Static method")

figure(2)
imagesc(X2_FE')
title("Dynamic method : decompositions")

figure(3)
imagesc(v')
title("Dynamic method : speed")