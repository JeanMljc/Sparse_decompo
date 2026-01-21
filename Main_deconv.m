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
SNR = 10;

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
kappa = 2;

% Build convolutive Gaussian dictionary
[A, mu, d_tail] = Dico_conv_gauss(N_lambda, std_problem, kappa);
M = length(mu);

%% Build ground truth signal & decomposition : Z_GT, X_GT 

% Ground-truth decompositions on the grid (does not handle oversampled grids)
X_GT = recover(amp_star, supp_star, M, N, d_tail);
Z_GT = A*X_GT;

% Error between Z_GT and Z_star
rMSEz_GT = rMSE(Z_GT,Z_star);

%% Solver hyperparameters CPLEX

% Param cplex
par_cplex.display = "off";
par_cplex.algo = 4;

%% Lp-norm approach: hyperparameters

par_LP.lambda = 10^-1; 
K_IRL1_LP = 10;

% Initialization of variables & weights
X0 = zeros(M,N);
W0 = ones(M,N);

%% Lp-norm approach : Call Primal-Dual CPLEX

X_LP = IRL1_CPLEX_LP(Z, A, W0, X0, par_LP, par_cplex, K_IRL1_LP, 0.9);
disp("End Lp-norm approach")

X_BPDN = X_LP{2};
X_ellp = X_LP{end};

Zre_LP = A*X_ellp;
rMSEx_LP = rMSE(X_ellp,X_GT);
rMSEz_LP = rMSE(Zre_LP,Z_star);

%% 1ST order : hyperparameters

par_1ST.omega = 10^-3; 
par_1ST.lambda = 10^-1; 
par_1ST.sdiag = 5*kappa;

K_IRWL1_OT1 = 10;

% Cost matrix
C_OT1 = cost_1st(M, 2);

% Initialization of variables
P0_1ST = zeros(M,M,N-1);
X0 = zeros(M,N);
W0 = ones(M,N);

%% 1ST order : Call Primal-Dual CPLEX

[X_1ST, P_1ST, w_IRL1, F1] = IRL1_CPLEX_1ST(Z, A, C_OT1, W0, X0, P0_1ST, par_1ST, par_cplex, K_IRWL1_OT1, 0.9);

disp("End 1ST ORDER")

X1_cvx = X_1ST{2};
X1_ncvx = X_1ST{end};

Zre_1ST = A*X1_ncvx;
rMSEx_1ST = rMSE(X1_ncvx,X_GT);
rMSEz_1ST = rMSE(Zre_1ST,Z_star);

%% 2ND order : hyperparameters

par_2ND.omega = 10^-2; 
par_2ND.lambda = 10^-1;
par_2ND.cmax = 2;

K_IRWL1_OT2 = 10;

% Cost tensor
C_OT2a = cost_2nd(M, 2); 

% Get indexes of the superband
idx_3D = find(cost_2nd(M,1) <= par_2ND.cmax);
par_2ND.ind = idx_3D;

% Superband cost tensor for OT term
C_2ND = C_OT2a(idx_3D);

% Initialization of algorithm
P0_2ND = [zeros((N-1)*M^2,1) ; repmat(zeros(size(C_2ND)),N-2,1)];

%% 2ND order : Call Primal-Dual CPLEX

[X_2ND, P_2ND] = IRL1_CPLEX_2ND(Z, A, C_2ND, W0, X0, P0_2ND, par_2ND, par_cplex, K_IRWL1_OT2, 0.9);

disp("End 2ND ORDER")

X2_cvx = X_2ND{2};
X2_ncvx = X_2ND{end};

Zre_2ND = A*X2_ncvx;
rMSEx_2ND = rMSE(X2_ncvx,X_GT);
rMSEz_2ND = rMSE(Zre_2ND,Z_star);

%% Plot part

figure(1)

subplot(221)
imagesc(X_GT',[0 3.5])
hold on
plot(kappa*(tra'-d_tail)-(kappa-1),1:N,'r','linewidth',1.5)
colorbar()
title("(a)")
axis equal
axis tight

subplot(222)
imagesc(X_ellp',[0 3.5])
hold on
plot(kappa*(tra'-d_tail)-(kappa-1),1:N,'r','linewidth',1.5)
colorbar()
title("(b)")
axis equal
axis tight

subplot(223)
imagesc(X1_ncvx',[0 3.5])
hold on
plot(kappa*(tra'-d_tail)-(kappa-1),1:N,'r','linewidth',1.5)
colorbar()
title("(c)")
axis equal
axis tight

subplot(224)
imagesc(X2_ncvx',[0 3.5])
hold on
plot(kappa*(tra'-d_tail)-(kappa-1),1:N,'r','linewidth',1.5)
colorbar()
title("(d)")
axis equal
axis tight

%% Plot figure 3
% 
% figure(1)
% 
% subplot(121)
% imagesc(X_GT',[0 3.5])
% hold on
% plot(ratio_super*(tra'-d_tail)-(ratio_super-1),1:N,'r','linewidth',1.5)
% colorbar()
% ylabel("Canal n")
% xlabel("Atom m")
% axis equal
% axis tight
% 
% subplot(122)
% imagesc(X_1ST{end}',[0 3.5])
% hold on
% plot(ratio_super*(tra'-d_tail)-(ratio_super-1),1:N,'r','linewidth',1.5)
% colorbar()
% xlabel("Atom m")
% axis equal
% axis tight
% 
% %%
% figure(2)
% 
% subplot(212)
% imagesc(X_2ND{end}',[0 3.5])
% hold on
% plot(ratio_super*(tra'-d_tail)-(ratio_super-1),1:N,'r','linewidth',1.5)
% colorbar()
% title("(d)")
% axis equal
% axis tight