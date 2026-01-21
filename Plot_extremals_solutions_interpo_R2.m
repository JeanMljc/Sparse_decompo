close all

% Use Interpo_R2.m to get successively X1 & X2 
% do not clear between

% Algo 1 & 2 give differents solutions X1 & X2
X3 = 0.5*(X1 + X2);
% Cost of X3 ? easy : the objective function is linear 

%% Plots

figure(1)
subplot(131)
imagesc(X1')
hold on 
plot(f,(1:N),'r','linewidth',2)
colorbar()
xlabel("Atome m")
ylabel("Canal n")
axis equal
axis tight

subplot(132)
imagesc(X2')
hold on 
plot(f,(1:N),'r','linewidth',2)
colorbar()
xlabel("Atome m")
ylabel("Canal n")
axis equal
axis tight

subplot(133)
imagesc(X3')
hold on 
plot(f,(1:N),'r','linewidth',2)
colorbar()
xlabel("Atome m")
ylabel("Canal n")
axis equal
axis tight
