clear
close all

% load path to subfolders
addpath(genpath('/home/jean/Documents/CODE_matlab/CODE_transmission/'))

%% Analytic expression of cost f

% function f(\alpha) contract by a factor K
f = @(x, K) (x <= 0.5/K) .* (2 * K * x) + (x > 0.5/K) .* (2*(1 - K * x));

%% Compress version of cost f : oversampled grid 

P = 2; % number of periode of f(\alpha) to plot 
K = 4; % \kappa goes from 1 to K
Nplot = 1000; % number of points to plot

fc = cell(1,K);
alphac = cell(1,K);

for k=(1:K)
    alpha = linspace(0, 1/k, Nplot);
    fc{k} = repmat(f(alpha,k) , 1, k*P);
    alphac{k} = linspace(0, P, length(fc{k}));
end 

%% Plot cost w.r.t alpha and theta

figure(1)
plot(alphac{1},fc{1})
hold on
plot(alphac{2},fc{2})
hold on
plot(alphac{3},fc{3})
ylabel("f(\alpha)")
xlabel("\alpha")
legend("\kappa=1","\kappa=2","\kappa=3")

% figure(2)
% plot(atan(alphac{1}),fc{1})
% hold on
% plot(atan(alphac{2}),fc{2})
% hold on
% plot(atan(alphac{3}),fc{3})
% ylabel("f(\theta)")
% xlabel("\theta = atan(\alpha)")
% legend("\kappa=1","\kappa=2","\kappa=3")
