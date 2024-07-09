%% compare model evidence
clc
clear
close all

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-4
    root_path = [root_path, str{i}, '/'];
end

addpath(genpath([root_path, 'tools/']))

% Loss Task
load([root_path, '/results/behav/loss_result_3e_al.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/loss_result_3o_al.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/loss_result_4_al.mat'], 'F_4_L1_params', 'F_4_L3_params')
F_3e_al_L1_params = F_3e_L1_params;
F_3o_al_L1_params = F_3o_L1_params;
F_4_al_L1_params = F_4_L1_params;
F_3e_al_L3_params = F_3e_L3_params;
F_3o_al_L3_params = F_3o_L3_params;
F_4_al_L3_params = F_4_L3_params;

load([root_path, '/results/behav/loss_result_3e_be.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/loss_result_3o_be.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/loss_result_4_be.mat'], 'F_4_L1_params', 'F_4_L3_params')
F_3e_be_L1_params = F_3e_L1_params;
F_3o_be_L1_params = F_3o_L1_params;
F_4_be_L1_params = F_4_L1_params;
F_3e_be_L3_params = F_3e_L3_params;
F_3o_be_L3_params = F_3o_L3_params;
F_4_be_L3_params = F_4_L3_params;

load([root_path, '/results/behav/loss_result_3e.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/loss_result_3o.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/loss_result_4.mat'], 'F_4_L1_params', 'F_4_L3_params')

[alpha,exp_r,xp,pxp,bor] = spm_BMS([F_3e_L1_params  F_3e_al_L1_params  F_3e_be_L1_params ...
    F_3o_L1_params  F_3o_al_L1_params  F_3o_be_L1_params ...
    F_4_L1_params  F_4_al_L1_params  F_4_be_L1_params]);
%{
[alpha,exp_r,xp,pxp,bor] = spm_BMS([F_3e_L1_params F_3e_L3_params F_3e_al_L1_params F_3e_al_L3_params F_3e_be_L1_params F_3e_be_L3_params ...
    F_3o_L1_params F_3o_L3_params F_3o_al_L1_params F_3o_al_L3_params F_3o_be_L1_params F_3o_be_L3_params ...
    F_4_L1_params F_4_L3_params F_4_al_L1_params F_4_al_L3_params F_4_be_L1_params F_4_be_L3_params]);
%}

    
figure; subplot(1,2,1); bar(pxp); ylabel('pxp'); 
title('Loss task');
%xticklabels({'3e+d', '3e+ad', '3e+al+d', '3e+al+ad', '3o+d', '3o+ad', '3o+al+d', '3o+al+ad', ... 
%    '4+d', '4+ad', '4+al+d', '4+al+ad'});
xlabel('Model')
pxp

% Reward task
load([root_path, '/results/behav/reward_result_3e_al.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/reward_result_3o_al.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/reward_result_4_al.mat'], 'F_4_L1_params', 'F_4_L3_params')
F_3e_al_L1_params = F_3e_L1_params;
F_3o_al_L1_params = F_3o_L1_params;
F_4_al_L1_params = F_4_L1_params;
F_3e_al_L3_params = F_3e_L3_params;
F_3o_al_L3_params = F_3o_L3_params;
F_4_al_L3_params = F_4_L3_params;

load([root_path, '/results/behav/reward_result_3e_be.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/reward_result_3o_be.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/reward_result_4_be.mat'], 'F_4_L1_params', 'F_4_L3_params')
F_3e_be_L1_params = F_3e_L1_params;
F_3o_be_L1_params = F_3o_L1_params;
F_4_be_L1_params = F_4_L1_params;
F_3e_be_L3_params = F_3e_L3_params;
F_3o_be_L3_params = F_3o_L3_params;
F_4_be_L3_params = F_4_L3_params;

load([root_path, '/results/behav/reward_result_3e.mat'], 'F_3e_L1_params', 'F_3e_L3_params')
load([root_path, '/results/behav/reward_result_3o.mat'], 'F_3o_L1_params', 'F_3o_L3_params')
load([root_path, '/results/behav/reward_result_4.mat'], 'F_4_L1_params', 'F_4_L3_params')

[alpha,exp_r,xp,pxp,bor] = spm_BMS([F_3e_L1_params  F_3e_al_L1_params  F_3e_be_L1_params  ...
    F_3o_L1_params  F_3o_al_L1_params  F_3o_be_L1_params  ...
    F_4_L1_params  F_4_al_L1_params  F_4_be_L1_params]);
%{
[alpha,exp_r,xp,pxp,bor] = spm_BMS([F_3e_L1_params F_3e_L3_params F_3e_al_L1_params F_3e_al_L3_params F_3e_be_L1_params F_3e_be_L3_params ...
    F_3o_L1_params F_3o_L3_params F_3o_al_L1_params F_3o_al_L3_params F_3o_be_L1_params F_3o_be_L3_params ...
    F_4_L1_params F_4_L3_params F_4_al_L1_params F_4_al_L3_params F_4_be_L1_params F_4_be_L3_params]);
%}
subplot(1,2,2); bar(pxp); ylabel('pxp'); 
title('Reward task');
%xticklabels({'3e+d', '3e+ad', '3e+al+d', '3e+al+ad', '3o+d', '3o+ad', '3o+al+d', '3o+al+ad', ... 
%    '4+d', '4+ad', '4+al+d', '4+al+ad'});
xlabel('Model')
pxp

set(gcf, 'color', 'w')

%% Comparing Estimated parameters from best model
% Loss
load('/home/dayoung/RL/ActInf/result/loss_est_2e.mat', 'Estimated_la_2e_L1', 'Estimated_eta_2e_L1')
figure; scatter(Estimated_la_2e_L1, Estimated_eta_2e_L1)
xlabel('la')
ylabel('eta')

% Reward
load('/home/dayoung/RL/ActInf/result/reward_est_4.mat', 'Estimated_alpha_4_L1', 'Estimated_rs_4_L1', 'Estimated_omega_4_L1', 'Estimated_eta_4_L1')
figure; scatter(Estimated_alpha_4_L1, Estimated_rs_4_L1)
xlabel('alpha')
ylabel('rs')

%% compare parameter estimation 
load('/home/dayoung/RL/ActInf/result/reward_est_4params.mat', 'Estimated_la_4', 'Estimated_la_4a', 'Estimated_rs_4', 'Estimated_rs_4a')
Estimated_la_4_re = Estimated_la_4;
Estimated_la_4a_re = Estimated_la_4a;
Estimated_rs_4_re = Estimated_rs_4;
Estimated_rs_4a_re = Estimated_rs_4a;

load('/home/dayoung/RL/ActInf/result/loss_est_4params.mat', 'Estimated_la_4', 'Estimated_la_4a', 'Estimated_rs_4', 'Estimated_rs_4a')
 
figure;
subplot(2,2,1)
scatter(Estimated_la_4_re, Estimated_la_4, 'filled')
lsline
title('Estimated la of inference model')
xlabel('reward task') 
ylabel('loss task') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4_re, Estimated_la_4);
text(3.5, 3, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3.5, 2.8, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,2)
scatter(Estimated_la_4a_re, Estimated_la_4a, 'filled')
lsline
title('Estimated la of learning model')
xlabel('reward task') 
ylabel('loss task') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4a_re, Estimated_la_4a);
text(3, 3.5, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3, 3.2, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,3)
scatter(Estimated_rs_4_re, Estimated_rs_4, 'filled')
lsline
title('Estimated rs of 4par')
xlabel('reward task') 
ylabel('loss task') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4_re, Estimated_rs_4);
text(4, 3.6, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(4, 3.5, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,4)
scatter(Estimated_rs_4a_re, Estimated_rs_4a, 'filled')
lsline
title('Estimated rs of 5par')
xlabel('reward task') 
ylabel('loss task') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4a_re, Estimated_rs_4a);
text(4.5, 5.4, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(4.5, 5.2, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

%%
figure; subplot(2,2,1)
scatter(Estimated_la_4 ./ Estimated_la_4_re, Estimated_rs_4_re ./ Estimated_rs_4, 'filled')
lsline
title('inference model')
xlabel('la loss / la reward') 
ylabel('rs reward / rs loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4 ./ Estimated_la_4_re, Estimated_rs_4_re ./ Estimated_rs_4);
text(1.4, 1.2, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(1.4, 1.1, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])

subplot(2,2,2)
scatter(Estimated_la_4a ./ Estimated_la_4a_re, Estimated_rs_4a_re ./ Estimated_rs_4a, 'filled')
lsline
title('learning model')
xlabel('la loss / la reward') 
ylabel('rs reward / rs loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4a ./ Estimated_la_4a_re, Estimated_rs_4a_re ./ Estimated_rs_4a);
text(3, 1.4, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3, 1.3, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])

subplot(2,2,3)
scatter(Estimated_la_4 ./ Estimated_la_4_re, Estimated_rs_4 ./ Estimated_rs_4_re, 'filled')
lsline
title('4par')
xlabel('la loss / la reward') 
ylabel('rs loss / rs reward') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4 ./ Estimated_la_4_re, Estimated_rs_4 ./ Estimated_rs_4_re);
text(1.3, 1, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(1.3, 0.95, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])

subplot(2,2,4)
scatter(Estimated_la_4a ./ Estimated_la_4a_re, Estimated_rs_4a ./ Estimated_rs_4a_re, 'filled')
lsline
title('5par')
xlabel('la loss / la reward') 
ylabel('rs loss / rs reward') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4a ./ Estimated_la_4a_re, Estimated_rs_4a ./ Estimated_rs_4a_re);
text(3, 1.6, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3, 1.4, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])
%%
figure; subplot(1,2,1)
scatter(Estimated_rs_4_re ./ Estimated_la_4_re, Estimated_rs_4 ./ Estimated_la_4, 'filled')
lsline
title('inference model')
xlabel('rs / la reward') 
ylabel('rs / la loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4_re ./ Estimated_la_4_re, Estimated_rs_4 ./ Estimated_la_4);
text(1.5, 1.6, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(1.5, 1.5, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])

subplot(1,2,2)
scatter(Estimated_rs_4a_re ./ Estimated_la_4a_re, Estimated_rs_4a ./ Estimated_la_4a, 'filled')
lsline
title('learning model')
xlabel('rs / la reward') 
ylabel('rs / la loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4a_re ./ Estimated_la_4a_re, Estimated_rs_4a ./ Estimated_la_4a);
text(3.4, 2.6, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3.4, 2.5, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])
%%
figure; subplot(1,2,1)
scatter(Estimated_rs_4_re ./ Estimated_la_4_re, Estimated_la_4 ./ Estimated_rs_4, 'filled')
lsline
title('inference model')
xlabel('rs / la reward') 
ylabel('la / rs loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4_re ./ Estimated_la_4_re, Estimated_la_4 ./ Estimated_rs_4);
text(1.5, 0.9, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(1.5, 0.85, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])

subplot(1,2,2)
scatter(Estimated_rs_4a_re ./ Estimated_la_4a_re, Estimated_la_4a ./ Estimated_rs_4a, 'filled')
lsline
title('learning model')
xlabel('rs / la reward') 
ylabel('la / rs loss') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_rs_4a_re ./ Estimated_la_4a_re, Estimated_la_4a ./ Estimated_rs_4a);
text(3.4, 0.8, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3.4, 0.75, ['p = ' num2str(round(Sig_alpha_2(1,2),5))])



%%
figure;
subplot(2,2,1)
scatter(Estimated_la_4_re, Estimated_rs_4_re, 'filled')
lsline
title('reward task inference model')
xlabel('loss aversive(la)') 
ylabel('reward sensitivity(rs)') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4_re, Estimated_rs_4_re);
text(3.5, 3.8, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3.5, 3.7, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,2)
scatter(Estimated_la_4a_re, Estimated_rs_4a_re, 'filled')
lsline
title('reward task learning model')
xlabel('loss aversive(la)') 
ylabel('reward sensitivity(rs)') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4a_re, Estimated_rs_4a_re);
text(3, 2.5, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3, 2.2, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,3)
scatter(Estimated_la_4, Estimated_rs_4, 'filled')
lsline
title('loss task inference model')
xlabel('loss aversive(la)') 
ylabel('reward sensitivity(rs)') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4, Estimated_rs_4);
text(3.5, 3.5, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3.5, 3.3, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])

subplot(2,2,4)
scatter(Estimated_la_4a, Estimated_rs_4a, 'filled')
lsline
title('loss task learning model')
xlabel('loss aversive(la)') 
ylabel('reward sensitivity(rs)') 
[Corr_alpha_2, Sig_alpha_2] = corrcoef(Estimated_la_4a, Estimated_rs_4a);
text(3, 3, ['r = ' num2str(round(Corr_alpha_2(1,2),4))])
text(3, 2.8, ['p = ' num2str(round(Sig_alpha_2(1,2),4))])
%% compare with survey data
load('/home/dayoung/RL/ActInf/result/reward_est_4params.mat')

Estimated_alpha_4a_re = Estimated_alpha_4a;
Estimated_beta_4a_re = Estimated_beta_4a;
Estimated_la_4a_re = Estimated_la_4a;
Estimated_rs_4a_re = Estimated_rs_4a;

load('/home/dayoung/RL/ActInf/result/reward_est_6params.mat')

Estimated_alpha_6_re = Estimated_alpha_6;
Estimated_beta_6_re = Estimated_beta_6;
Estimated_la_6_re = Estimated_la_6;
Estimated_rs_6_re = Estimated_rs_6;

corr(Estimated_alpha_4a_re, Estimated_alpha_6_re)
corr(Estimated_beta_4a_re, Estimated_beta_6_re)
corr(Estimated_la_4a_re, Estimated_la_6_re)
corr(Estimated_rs_4a_re, Estimated_rs_6_re)

%%
load('/home/dayoung/RL/ActInf/result/loss_est_4params.mat')

%%
loss_5 = Estimated_la_3./Estimated_rs_3;
loss_4 = Estimated_la_2./Estimated_rs_2;
reward_5 = Estimated_la_3_re./Estimated_rs_3_re;
reward_4 = Estimated_la_2_re./Estimated_rs_2_re;
%
[Corr1, p_val1] = corrcoef(loss_5, phq_sum');
[Corr2, p_val2] = corrcoef(loss_4, phq_sum');
[Corr3, p_val3] = corrcoef(reward_5, phq_sum');
[Corr4, p_val4] = corrcoef(reward_4, phq_sum');

for kk = 1:size(phq,1)
    kk
    [Corr1, p_val1] = corrcoef(loss_5, phq(kk,:)')
    [Corr2, p_val2] = corrcoef(loss_4, phq(kk,:)')
    [Corr3, p_val3] = corrcoef(reward_5, phq(kk,:)')
    [Corr4, p_val4] = corrcoef(reward_4, phq(kk,:)')
end
%%
%
loss_5 = Estimated_la_3./Estimated_rs_3_re;
loss_4 = Estimated_la_2./Estimated_rs_2_re;
 
%
for kk = 1:size(phq,1)
    kk
    [Corr1, p_val1] = corrcoef(loss_5, phq(kk,:)')
    [Corr2, p_val2] = corrcoef(loss_4, phq(kk,:)')
end