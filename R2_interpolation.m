clear
close all

% load path to subfolders
addpath(genpath('/home/jean/Documents/CODE_matlab/CODE_transmission/'))

% load path to CPLEX solver 
addpath('/home/jean/CPLEX_1210/cplex/matlab/x86-64_linux')

%% Build data to interpolate

N = 11; % number of channels to interpolate 
alpha  = 0.7; % slop of discrete line
theta = atan(alpha); % angle associated to alpha

y_1 = 5;
y_N = round(y_1 + (N-1)*alpha);

% Value of oversampling factor 
kappa = [1 2 3];

%% Build two marginals to interpolate

Marg = cell(1,length(kappa)); % decomposition to interpolate 
f = cell(1,length(kappa)); % continuous line  
C = cell(1,length(kappa)); % OT cost
X0 = cell(1,length(kappa)); % Initialization of algo 

for k=1:length(kappa)
    M = kappa(k)*(y_N + 5); % size of oversampled dictionary 
    Marg{k} = full(sparse([kappa(k)*y_1 kappa(k)*y_N],[1 2],1, M, 2));
    f{k} = linspace(kappa(k)*y_1,kappa(k)*y_N,N);
    C{k} = cost_2nd(M,2); %squared cost
    X0{k} = zeros(M,N);
end

%% Solver hyperparameters CPLEX

% Param cplex
par_cplex.display = "off";
par_cplex.algo = 0;

% Call solver CPLEX to solve

X = cell(1,length(kappa)); % solution
F_min = cell(1,length(kappa)); % objective value

% CPLEX interpolation 

for k=1:length(kappa)
    [~, ~, X{k}, F_min{k}, ~, ~] = Algo_interpo_R2(C{k}, Marg{k}, X0{k}, N, par_cplex);
    disp(k)
end

%% Discrete line on the grid 

% For kappa=1 only 

line = floor((1:N)*alpha) + y_1; % line on grid 
X_line = full(sparse((1:N)', line, ones(N,1), N, (y_N + 5)))';

%% Plots

figure(1)

subplot(131)
imagesc(X{1}')
hold on 
plot(f{1},(1:N),'r','linewidth',1)
colorbar()
title("\kappa = 1")
xlabel("Atome m")
ylabel("Canal n")

subplot(132)
imagesc(X{2}')
hold on 
plot(f{2},(1:N),'r','linewidth',1)
colorbar()
title("\kappa = 2")
xlabel("Atome m")
ylabel("Canal n")

subplot(133)
imagesc(X{3}')
hold on 
plot(f{3},(1:N),'r','linewidth',1)
colorbar()
title("\kappa = 3")
xlabel("Atome m")
ylabel("Canal n")

% figure(2)
% 
% subplot(121)
% imagesc(X{1}')
% hold on 
% plot(f{1},(1:N),'r','linewidth',1.5)
% colorbar()
% title("Interpolation solution")
% xlabel("Atome m")
% ylabel("Canal n")
% 
% subplot(122)
% imagesc(X_line')
% hold on 
% plot(f{1},(1:N),'r','linewidth',1.5)
% colorbar()
% title("Discrete line")
% xlabel("Atome m")
% ylabel("Canal n")


