close all

% This script aim to plot transport displacement from OT plans

% First execute Main_deconv_ncvx_cross_traj.m 
% to get X_1ST and X_2ND

%% 1ST order approach 

% Get sequence of decompositions
Xp_1ST = X_1ST{end};

% Get mass displacement
sub_OT_1ST = Get_OT_mass_disp(P_1ST,0.8);
ART_1ST = Get_arrow(sub_OT_1ST);

figure(1)
imagesc(Xp_1ST',[0 4])
colorbar()
hold on
q = quiver(ART_1ST(:,2),ART_1ST(:,1),ART_1ST(:,4),ART_1ST(:,3),0,'r','LineWidth',2);
q.MaxHeadSize = 0.1;
title("R1")
axis equal
axis tight
xlabel("m atoms")
ylabel("n channels")

%% 2ND order approach 

Xp_2ND = X_2ND{end};

sub_OT_2ND = Get_OT_mass_disp(P_2ND, 0.8);
ART_2ND = Get_arrow(sub_OT_2ND);

figure(2)
imagesc(Xp_2ND',[0 4])
colorbar()
hold on
q = quiver(ART_2ND(:,2),ART_2ND(:,1),ART_2ND(:,4),ART_2ND(:,3),0,'r','LineWidth',2);
q.MaxHeadSize = 0.1;
title("R2")
axis equal
axis tight
xlabel("m atoms")
ylabel("n channels")