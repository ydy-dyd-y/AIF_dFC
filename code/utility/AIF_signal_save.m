clear
clc
close all

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

addpath(genpath([root_path, 'tools/']))
%%
tp = 231;
TR = 1.5;

load([root_path, 'data/behav/action_fix_time.mat'], 'action_time', 'fix_time');
action_time(isnan(action_time)) = 1.5;
load([root_path, 'data/behav/behav_measure.mat'], 'lo_feed', 'lo_cue', 're_feed', 're_cue');

feed_time = [lo_feed, re_feed] - 7.5;
cue_time = [lo_cue, re_cue] - 7.5;

load([root_path, 'data/behav/feedback_win_or_lose.mat'], 'feedback_win_or_lose')
feedback_win_or_lose(:,[20, 32]) = [];

load([root_path, 'results/behav/loss_result_3o_al.mat'], 'GMDP_3o_L1')
GMDP_loss = GMDP_3o_L1;
load([root_path, 'results/behav/reward_result_4_al.mat'], 'GMDP_4_L1')
GMDP_reward = GMDP_4_L1;

rm_fields = {'qln', 'errn', 'ERRn'};
for si = 1:numel(GMDP_loss)
    GMDP_loss{si} = rmfield(GMDP_loss{si},rm_fields);
    GMDP{si, 1} = [GMDP_loss{si}, GMDP_reward{si}];
end

num_trial = size(action_time, 1);
Nsubj = size(action_time, 2);

x = round(linspace(0, 348, 3480),1);
v = hrf(x);

rm_sub1 = [20, 32];
rm_sub2 = [22, 36, 44, 57];
action_time(:, rm_sub1) = []; action_time(:, rm_sub2) = [];
fix_time(:, rm_sub1) = []; fix_time(:, rm_sub2) = [];

feedback_win_or_lose(:, rm_sub2) = [];
%% Precision
Nsubj_neuro = Nsubj - 2;

GMDP_dn = cell(Nsubj_neuro,num_trial);
GMDP_wn = cell(Nsubj_neuro,num_trial);
ii = 1;
for si = 1:Nsubj  % 60: different with number of imaging data(58)
    if (si ~= 20) && (si ~= 32)
        for ti = 1:num_trial
            if ti < 41
                GMDP_dn{ii,ti} = GMDP_loss{si,1}(ti).dn;   % [iter(16) * T(4) x 1]
                GMDP_wn{ii,ti} = GMDP_loss{si,1}(ti).wn; 
            else
                GMDP_dn{ii,ti} = GMDP_reward{si,1}(ti - 40).dn;   % [iter(16) * T(4) x 1]
                GMDP_wn{ii,ti} = GMDP_reward{si,1}(ti - 40).wn; 
            end
        end
        ii = ii + 1;
    else
        continue
    end
end

W_dn = zeros(length(GMDP_dn{1,1}),num_trial,Nsubj_neuro);  % [iter(=16) * T(=2) x 40 x 58]
W_wn = zeros(length(GMDP_wn{1,1}),num_trial,Nsubj_neuro);  % [iter(=16) * T(=2) x 40 x 58]
for si = 1:Nsubj_neuro
    for ti =1:num_trial
        w_dn(:,ti) = GMDP_dn{si,ti};
        w_wn(:,ti) = GMDP_wn{si,ti};
    end
    W_dn(:,:,si) = w_dn;
    W_wn(:,:,si) = w_wn;
end
W_dn(:, :, rm_sub2) = [];
W_wn(:, :, rm_sub2) = [];

Precision = zeros(size(W_dn, 3), tp*TR*10);
Precision_hrf = [];
max_dn = [];

for si = 1:size(Precision, 1)
    max_dn_indi = [];
    for ti = 1:num_trial
        if ti == 1 || ti == 41
            max_dn_indi = [max_dn_indi, 0.125];
        else
            max_dn_indi = [max_dn_indi, max(W_dn(:, ti - 1, si))];
            for wti = 1:size(W_dn,1)
                if (wti > 16) && (W_dn(wti, ti - 1, si)>=0)
                    Precision(si, round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10):round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10)+9) = repelem(W_dn(wti,ti - 1,si),10);
                end
            end
        end
    end
    max_dn = [max_dn; max_dn_indi];
    Precision_hrf_tmp = conv(Precision(si,:),v);
    Precision_hrf(si, :) = Precision_hrf_tmp;        
end
Precision_hrf = Precision_hrf(:, 1:3465);
%%
Precision2 = zeros(size(W_wn, 3), tp*TR*10);
Precision_hrf2 = [];
max_wn = [];
for si = 1:size(Precision2, 1)
    max_wn_indi = [];
    for ti = 1:length(feed_time)
        if ti == 1 || ti == 41
            for wti = 1:16
                Precision2(si, round((feed_time(ti)+(1.0/16)*(wti-17))*10):round((feed_time(ti)+(1.0/16)*(wti-17))*10)+9) = repelem(W_wn(wti,ti,si),10);
            end
            max_wn_indi = [max_wn_indi, 1.25];
        else
            max_wn_indi = [max_wn_indi, max(W_wn(end, ti - 1, si))];
            for wti = 1:size(W_wn,1)
                if (wti < 17)  
                    Precision2(si, round((feed_time(ti)+(1.0/16)*(wti-17))*10):round((feed_time(ti)+(1.0/16)*(wti-17))*10)+9) = repelem(W_wn(wti,ti,si),10);                    
                else
                    Precision2(si, round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10):round((cue_time(ti)+(action_time(ti,si)/16)*(wti-1))*10)+9) = repelem(W_wn(wti,ti - 1,si),10);
                end
            end
        end
    end
    max_wn = [max_wn; max_wn_indi];
    Precision_hrf_tmp = conv(Precision2(si,:),v);
    Precision_hrf2(si, :) = Precision_hrf_tmp;        
end
Precision_hrf2 = Precision_hrf2(:, 1:3465);
%%
sub_id = 17;
%figure
subplot(2,1,1), area(1:length(Precision(sub_id,:)), Precision(sub_id,:))
subplot(2,1,2), plot(Precision_hrf(sub_id,:)) %area(1:length(Precision2(sub_id,:)), Precision2(sub_id,:))
set(gcf, 'color', 'w')

%% precision and action time correlation
corr_loss_wn = []; corr_reward_wn = [];
for si = 1:size(max_wn, 1)
    tmp = corrcoef(max_wn(si, 1:40), 1./action_time(1:40,si));
    corr_loss_wn = [corr_loss_wn, tmp(1,2)];

    tmp = corrcoef(max_wn(si, 41:80), 1./action_time(41:80,si));
    corr_reward_wn = [corr_reward_wn, tmp(1,2)];
end

figure; subplot(2,2,1); histogram(corr_loss_wn)
xline(nanmean(corr_loss_wn), 'lineWidth', 1.5); ylim([0, 20])
subplot(2,2,2); histogram(corr_reward_wn)
xline(nanmean(corr_reward_wn), 'lineWidth', 1.5); ylim([0, 20])
nanmean(corr_loss_wn)
nanmean(corr_reward_wn)

sub_id = 51; 
subplot(2,2,3); scatter(max_wn(sub_id, 1:40), action_time(1:40,sub_id))
subplot(2,2,4); scatter(max_wn(sub_id, 41:80), action_time(41:80,sub_id))
set(gcf, 'color', 'w')

%% precision gradient and action time correlation
corr_loss_dn = []; corr_reward_dn = [];
for si = 1:size(max_dn, 1)
    tmp = corrcoef(max_dn(si, 1:40), 1./action_time(1:40,si));
    corr_loss_dn = [corr_loss_dn, tmp(1,2)];

    tmp = corrcoef(max_dn(si, 41:80), 1./action_time(41:80,si));
    corr_reward_dn = [corr_reward_dn, tmp(1,2)];
end

figure; subplot(2,2,1); histogram(corr_loss_dn, 15)
xline(nanmean(corr_loss_dn), 'lineWidth', 1.5); ylim([0, 12])
subplot(2,2,2); histogram(corr_reward_dn, 15)
xline(nanmean(corr_reward_dn), 'lineWidth', 1.5); ylim([0, 12])
nanmean(corr_loss_dn)
nanmean(corr_reward_dn)

sub_id = 51;
subplot(2,2,3); scatter(max_dn(sub_id, 1:40), action_time(1:40,sub_id))
subplot(2,2,4); scatter(max_dn(sub_id, 41:80), action_time(41:80,sub_id))
set(gcf, 'color', 'w')

%{
action_time_tmp = action_time';
corrcoef(reshape(max_dn, [1, numel(max_dn)]), reshape(action_time_tmp, [1, numel(action_time_tmp)]))
figure; scatter(reshape(max_dn, [1, numel(max_dn)]), reshape(action_time_tmp, [1, numel(action_time_tmp)]))
set(gcf, 'color', 'w')
%}
%% correlation of precision policy and precision of prediction
load([root_path, 'results/behav/hgf.mat'], 'est_hgf_lore')
varL1 = []; varR1 = [];
varL2 = []; varR2 = [];
varL3 = []; varR3 = [];

for si = 1:Nsubj
    if (si ~= 20) && (si ~= 32)
        % variance of predictions
        varL1 = [varL1, est_hgf_lore{si,1}.traj.sahat(:,1)];
        varR1 = [varR1, est_hgf_lore{si,2}.traj.sahat(:,1)];
        varL2 = [varL2, est_hgf_lore{si,1}.traj.sahat(:,2)];
        varR2 = [varR2, est_hgf_lore{si,2}.traj.sahat(:,2)];
        varL3 = [varL3, est_hgf_lore{si,1}.traj.sahat(:,3)];
        varR3 = [varR3, est_hgf_lore{si,2}.traj.sahat(:,3)];
    end
end
varL1(:,rm_sub2) = []; varR1(:,rm_sub2) = [];
varL2(:,rm_sub2) = []; varR2(:,rm_sub2) = [];
varL3(:,rm_sub2) = []; varR3(:,rm_sub2) = [];

corr_loss_wn = []; corr_reward_wn = [];
corr_loss_dn = []; corr_reward_dn = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    tmp = corrcoef(max_wn(si, 1:40), 1./varL1(:,si));
    corr_loss_wn = [corr_loss_wn, tmp(1,2)];
    tmp = corrcoef(max_wn(si, 41:80), 1./varR1(:,si));
    corr_reward_wn = [corr_reward_wn, tmp(1,2)];

    tmp = corrcoef(max_dn(si, 1:40), 1./varL1(:,si));
    corr_loss_dn = [corr_loss_dn, tmp(1,2)];
    tmp = corrcoef(max_dn(si, 41:80), 1./varR1(:,si));
    corr_reward_dn = [corr_reward_dn, tmp(1,2)];
end

figure; subplot(2,2,1); histogram(corr_loss_wn, 15)
hold on; xline(nanmean(corr_loss_wn), 'lineWidth', 1.5); ylim([0, 12])
subplot(2,2,2); histogram(corr_reward_wn, 15)
hold on; xline(nanmean(corr_reward_wn), 'lineWidth', 1.5); ylim([0, 12])
nanmean(corr_loss_wn)
nanmean(corr_reward_wn)

sub_id = 7;
subplot(2,2,3); scatter(max_wn(sub_id, 1:40), 1./varL1(:,sub_id))
subplot(2,2,4); scatter(max_wn(sub_id, 41:80), 1./varR1(:,sub_id))
sgtitle('Precision (gamma)')
set(gcf, 'color', 'w')

figure; subplot(2,2,1); histogram(corr_loss_dn, 15)
hold on; xline(nanmean(corr_loss_dn), 'lineWidth', 1.5); ylim([0, 12])
subplot(2,2,2); histogram(corr_reward_dn, 15)
hold on; xline(nanmean(corr_reward_dn), 'lineWidth', 1.5); ylim([0, 12])
nanmean(corr_loss_dn)
nanmean(corr_reward_dn)

sub_id = 7;
subplot(2,2,3); scatter(max_dn(sub_id, 1:40), 1./varL1(:,sub_id))
subplot(2,2,4); scatter(max_dn(sub_id, 41:80), 1./varR1(:,sub_id))
sgtitle('rate of Precision (phasic dopa)')
set(gcf, 'color', 'w')

%%
load([root_path, 'data/behav/PL_idx.mat'])
prec_wrongL16 = max_wn(wrongL16, 16);
prec_wrongL24 = max_wn(wrongL24, 24);
prec_wrongR23 = max_wn(wrongR23, 23 + 40);
prec_wrongR37 = max_wn(wrongR37, 37 + 40);

prec_rightL16 = max_wn(rightL16, 16);
prec_rightL24 = max_wn(rightL24, 24);
prec_rightR23 = max_wn(rightR23, 23 + 40);
prec_rightR37 = max_wn(rightR37, 37 + 40);

mean(prec_wrongL16)
mean(prec_rightL16)
mean(prec_wrongL24)
mean(prec_rightL24)
mean(prec_wrongR23)
mean(prec_rightR23)
mean(prec_wrongR37)
mean(prec_rightR37)

mean(prec_wrongL16(find(prec_wrongL16 < 1/1.5)))
mean(prec_rightL16(find(prec_rightL16 < 1/1.5)))
mean(prec_wrongL24(find(prec_wrongL24 < 1/1.5)))
mean(prec_rightL24(find(prec_rightL24 < 1/1.5)))
mean(prec_wrongR23(find(prec_wrongR23 < 1/1.5)))
mean(prec_rightR23(find(prec_rightR23 < 1/1.5)))
mean(prec_wrongR37(find(prec_wrongR37 < 1/1.5)))
mean(prec_rightR37(find(prec_rightR37 < 1/1.5)))

figure;
subplot(2, 2, 1); histogram(prec_wrongL16, 6); hold on; histogram(prec_rightL16, 6)
subplot(2, 2, 2); histogram(prec_wrongL24, 6); hold on; histogram(prec_rightL24, 6)
subplot(2, 2, 3); histogram(prec_wrongR23, 6); hold on; histogram(prec_rightR23, 6)
subplot(2, 2, 4); histogram(prec_wrongR37, 6); hold on; histogram(prec_rightR37, 6)

%% prediction error  
GMDP_vn = cell(Nsubj_neuro,num_trial);
ii = 1;
for si = 1:numel(GMDP)  % 60: different with number of imaging data(58)
    if (si ~= 20) && (si ~= 32)
        for ti = 1:num_trial
            GMDP_vn{ii,ti} = GMDP{si,1}(ti).vn{1,1}(:,1,2,2);   % (iteration, 1=fig.A better state, tau=2, t=2) -> [iter(16) x 1]
        end
        ii = ii + 1;
    else
        continue
    end
end

V = zeros(length(GMDP_vn{1,1}),num_trial,Nsubj_neuro);  % [iter(=16) x 40 x 58]

for si = 1:Nsubj_neuro
    for ti =1:num_trial
        v_tmp(:,ti) = GMDP_vn{si,ti};  % v_tmp : [iter(=16) x 40]
    end
    V(:,:,si) = v_tmp;
end

prediction_error = zeros(Nsubj_neuro, tp*TR*10);
for si = 1:Nsubj_neuro
    for ti = 1:length(feed_time)
        for vti = 1:size(V,1)
            prediction_error(si, round((feed_time(ti)+(1.0/16)*vti)*10):round((feed_time(ti)+(1.0/16)*vti)*10)+9) = repelem(V(vti,ti,si),10);
        end
    end
end

Pos_prediction_error = prediction_error;
Neg_prediction_error = prediction_error;
Signed_RPE = prediction_error;

Pos_prediction_error(find(Pos_prediction_error < 0)) = 0;   % PPE
Neg_prediction_error(find(Neg_prediction_error > 0)) = 0;   % NPE
Neg_prediction_error = (Neg_prediction_error) * (-1);
UnSigned_RPE = Pos_prediction_error + Neg_prediction_error;

for si = 1:Nsubj_neuro
    Pos_prediction_error_hrf_tmp = conv(Pos_prediction_error(si,:),v);
    Pos_prediction_error_hrf(si, :) = Pos_prediction_error_hrf_tmp;
    Neg_prediction_error_hrf_tmp = conv(Neg_prediction_error(si,:),v);
    Neg_prediction_error_hrf(si, :) = Neg_prediction_error_hrf_tmp;
    
    Signed_RPE_hrf_tmp = conv(Signed_RPE(si,:),v);
    Signed_RPE_hrf(si, :) = Signed_RPE_hrf_tmp;
    UnSigned_RPE_hrf_tmp = conv(UnSigned_RPE(si,:),v);
    UnSigned_RPE_hrf(si, :) = UnSigned_RPE_hrf_tmp;
end
Pos_prediction_error_hrf = Pos_prediction_error_hrf(:, 1:3465);
Neg_prediction_error_hrf = Neg_prediction_error_hrf(:, 1:3465);
Signed_RPE_hrf = Signed_RPE_hrf(:, 1:3465);
UnSigned_RPE_hrf = UnSigned_RPE_hrf(:, 1:3465);

% 58 -> 54
Pos_prediction_error_hrf(rm_sub2, :) = [];   % [Nsubj x 3414] -> [Nsubj x (tp*TR*10)]
Pos_prediction_error(rm_sub2, :) = []; 
Neg_prediction_error_hrf(rm_sub2, :) = [];   
Neg_prediction_error(rm_sub2, :) = []; 

Signed_RPE_hrf(rm_sub2, :) = [];   
Signed_RPE(rm_sub2, :) = []; 
UnSigned_RPE_hrf(rm_sub2, :) = [];   
UnSigned_RPE(rm_sub2, :) = []; 
%%
save([root_path, 'data/behav/signal_data.mat'], 'Precision_hrf', 'Precision', ...
    'Pos_prediction_error_hrf', 'Pos_prediction_error', 'Neg_prediction_error_hrf', 'Neg_prediction_error', ...
    'Signed_RPE_hrf', 'Signed_RPE', 'UnSigned_RPE_hrf', 'UnSigned_RPE')

%% correlation between PPE&NPE and precision of prediction
PPE_max = zeros(Nsubj_neuro, 80); 
NPE_max = zeros(Nsubj_neuro, 80);
for si = 1:Nsubj_neuro
    for ti = 1:80
        if max(squeeze(V(:, ti, si))) > 0
            PPE_max(si, ti) = max(squeeze(V(:, ti, si)));
        else
            NPE_max(si, ti) = max(squeeze(abs(V(:, ti, si))));
        end
    end
end
PPE_max(rm_sub2, :) = [];
NPE_max(rm_sub2, :) = [];

corr_loss_PPE = []; corr_reward_PPE = [];
corr_loss_NPE = []; corr_reward_NPE = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    tmp = corrcoef(PPE_max(si, 1:40), 1./varL2(:,si));
    corr_loss_PPE = [corr_loss_PPE, tmp(1,2)];
    tmp = corrcoef(PPE_max(si, 41:80), 1./varR2(:,si));
    corr_reward_PPE = [corr_reward_PPE, tmp(1,2)];

    tmp = corrcoef(NPE_max(si, 1:40), 1./varL2(:,si));
    corr_loss_NPE = [corr_loss_NPE, tmp(1,2)];
    tmp = corrcoef(NPE_max(si, 41:80), 1./varR2(:,si));
    corr_reward_NPE = [corr_reward_NPE, tmp(1,2)];
end

figure; subplot(2,2,1); histogram(corr_loss_PPE, 15)
hold on; xline(nanmean(corr_loss_PPE), 'lineWidth', 1.5); ylim([0, 12])
subplot(2,2,2); histogram(corr_reward_PPE, 15)
hold on; xline(nanmean(corr_reward_PPE), 'lineWidth', 1.5); ylim([0, 12])
nanmean(corr_loss_PPE)
nanmean(corr_reward_PPE)

sub_id = 7;
subplot(2,2,3); scatter(max_wn(sub_id, 1:40), 1./varL1(:,sub_id))
subplot(2,2,4); scatter(max_wn(sub_id, 41:80), 1./varR1(:,sub_id))
sgtitle('Precision (gamma)')
set(gcf, 'color', 'w')

figure; subplot(2,2,1); histogram(corr_loss_dn, 15)
hold on; xline(nanmean(corr_loss_dn), 'lineWidth', 1.5); ylim([0, 12])
subplot(2,2,2); histogram(corr_reward_dn, 15)
hold on; xline(nanmean(corr_reward_dn), 'lineWidth', 1.5); ylim([0, 12])
nanmean(corr_loss_dn)
nanmean(corr_reward_dn)

sub_id = 7;
subplot(2,2,3); scatter(max_dn(sub_id, 1:40), 1./varL1(:,sub_id))
subplot(2,2,4); scatter(max_dn(sub_id, 41:80), 1./varR1(:,sub_id))
sgtitle('rate of Precision (phasic dopa)')
set(gcf, 'color', 'w')
%%
for sub_id = 1:Nsubj_neuro - length(rm_sub2)
    if rem(sub_id, 10) == 1
        figure
    end
    if rem(sub_id, 10) == 0
        subplot(5,2,10), area(1:length(Pos_prediction_error(sub_id,:)), Pos_prediction_error(sub_id,:))
    else
        subplot(5,2,rem(sub_id, 10)), area(1:length(Pos_prediction_error(sub_id,:)), Pos_prediction_error(sub_id,:))
    end
end

%% prior state belief (d)
GMDP_d = cell(Nsubj_neuro,num_trial);
ii = 1;
for si = 1:numel(GMDP)  % 60: different with number of imaging data(58)
    if (si ~= 20) && (si ~= 32)
        for ti = 1:num_trial
            GMDP_d{ii,ti} = GMDP{si,1}(ti).d{1,1};  % 1 = fig.A better state
        end
        ii = ii + 1;
    else
        continue
    end
end

d_belief = zeros(length(GMDP_d{1,1}),num_trial,Nsubj_neuro);  % [2 x 40 x 58]

for si = 1:Nsubj_neuro
    for ti =1:num_trial
        d_belief_tmp(:,ti) = GMDP_d{si,ti};  % d_belief_tmp : [2 x 40]
    end
    d_belief(:,:,si) = d_belief_tmp;
end
d_belief(:,:,rm_sub2) =[];

for si = 1:10 %Nsubj_neuro-length(rm_sub2)
    if rem(si, 10) == 1
        figure;
    end
    if rem(si, 10) == 0
        subplot(10, 1, 10)
        plot(d_belief(1,:,si))
        hold on
        plot(d_belief(2,:,si))
        hold on
        %plot(d_belief(1,:,si)./d_belief(2,:,si), 'k')
    else
        subplot(10, 1, rem(si, 10))
        plot(d_belief(1,:,si))
        hold on
        plot(d_belief(2,:,si))
        hold on
        %plot(d_belief(1,:,si)./d_belief(2,:,si), 'k')   % corresponding to precision of state prediction
    end
end

%% state belief (xn)
GMDP_xn2 = cell(Nsubj_neuro,num_trial);
GMDP_xn1 = cell(Nsubj_neuro,num_trial);
ii = 1;
for si = 1:numel(GMDP)  % 60: different with number of imaging data(58)
    if (si ~= 20) && (si ~= 32)
        for ti = 1:num_trial
            GMDP_xn1{ii,ti} = GMDP{si,1}(ti).xn{1,1}(:,1,2,1); 
            GMDP_xn2{ii,ti} = GMDP{si,1}(ti).xn{1,1}(:,1,2,2);   % (iteration, 1=fig.A better state, tau=2, t=2) -> [iter(16) x 1]
        end
        ii = ii + 1;
    else
        continue
    end
end

S1 = zeros(length(GMDP_xn1{1,1}),num_trial,Nsubj_neuro);  % [iter(=16) x 40 x 58]
S2 = zeros(length(GMDP_xn2{1,1}),num_trial,Nsubj_neuro);  % [iter(=16) x 40 x 58]

for si = 1:Nsubj_neuro
    for ti =1:num_trial
        s1_tmp(:,ti) = GMDP_xn1{si,ti};  % v_tmp : [iter(=16) x 40]
        s2_tmp(:,ti) = GMDP_xn2{si,ti};  % v_tmp : [iter(=16) x 40]
    end
    S1(:,:,si) = s1_tmp;
    S2(:,:,si) = s2_tmp;
end

state_prediction = zeros(Nsubj_neuro, tp*TR*10);
for si = 1:Nsubj_neuro
    for ti = 1:length(feed_time)
        for vti = 1:size(S1,1)
            state_prediction(si, round((cue_time(ti)+(1.0/16)*vti)*10):round((cue_time(ti)+(1.0/16)*vti)*10)+9) = repelem(S1(vti,ti,si),10);
            state_prediction(si, round((feed_time(ti)+(1.0/16)*vti)*10):round((feed_time(ti)+(1.0/16)*vti)*10)+9) = repelem(S2(vti,ti,si),10);
        end
    end
end

for si = 1:Nsubj_neuro
    state_prediction_hrf(si, :) = conv(state_prediction(si,:),v);
end
state_prediction_hrf = state_prediction_hrf(:, 1:3465);

% 58 -> 54
rm_sub2 = [22, 36, 44, 57];
state_prediction(rm_sub2, :) = [];   % [Nsubj x 3414] -> [Nsubj x (tp*TR*10)]
state_prediction_hrf(rm_sub2, :) = []; 

sub_id = 19;
figure
subplot(1,2,1), area(1:length(state_prediction(sub_id,:)), state_prediction(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('State_prediction');t.FontSize = 8;

subplot(1,2,2), plot(1:length(state_prediction_hrf(sub_id,:)), state_prediction_hrf(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('State_prediction hrf');t.FontSize = 8; 

save([root_path, 'data/behav/signal_data.mat'], 'state_prediction', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'state_prediction_hrf', '-append')

%% starting task(ST)
ST = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    %for ti = [1 41]
        %ST(si, round((cue_time(1) - fix_time(1, si))*10):round((cue_time(1) - fix_time(1, si))*10) + 9) = 1;
        ST(si, round((cue_time(1) - fix_time(1, si))*10):round(cue_time(1))*10 + 9) = 1;
        % ST(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1;  % maybe cue is not the start of the task
        %ST(si, round((19-7.5)*10):round((19-7.5)*10) + 9) = 1;
        %ST(si, round((179-7.5)*10):round((179-7.5)*10) + 9) = 1;
        %ST(si, round((cue_time(41) - fix_time(41, si))*10):round((cue_time(41) - fix_time(41, si))*10) + 9) = 1;
        ST(si, round((cue_time(41) - fix_time(41, si))*10):round(cue_time(41))*10 + 9) = 1;
    %end
end

for si = 1:Nsubj_neuro - length(rm_sub2)
    ST_hrf_tmp = conv(ST(si,:),v);
    ST_hrf(si, :) = ST_hrf_tmp;
end
ST_hrf = ST_hrf(:, 1:3465);

%% cue visual stimulus(CVS)
CVS = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(cue_time)
        CVS(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1;
    end
end

for si = 1:Nsubj_neuro - length(rm_sub2)
    CVS_hrf_tmp = conv(CVS(si,:),v);
    CVS_hrf(si, :) = CVS_hrf_tmp;
end
CVS_hrf = CVS_hrf(:, 1:3465);

%% feedback visulat stimulus(FVS_w and FVS_l)
FVS_w = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_l = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);

for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(feed_time)
        if feedback_win_or_lose(ti, si) == 1
            FVS_w(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; %repelem(1,10);
        else
            FVS_l(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; %repelem(1,10);
        end
    end
end

FVS = FVS_w + FVS_l;
for si = 1:Nsubj_neuro - length(rm_sub2)
    FVS_w_hrf_tmp = conv(FVS_w(si,:),v);
    FVS_w_hrf(si, :) = FVS_w_hrf_tmp;
    FVS_l_hrf_tmp = conv(FVS_l(si,:),v);
    FVS_l_hrf(si, :) = FVS_l_hrf_tmp;    
    FVS_hrf_tmp = conv(FVS(si,:),v);
    FVS_hrf(si, :) = FVS_hrf_tmp;     
end
FVS_w_hrf = FVS_w_hrf(:, 1:3465);
FVS_l_hrf = FVS_l_hrf(:, 1:3465);
FVS_hrf = FVS_hrf(:, 1:3465);

%%
VS = CVS + FVS;
VS(VS == 2) = 1;

for si = 1:size(VS,1)
    VS1_hrf_tmp = conv(VS(si,:),v);
    VS1_hrf(si, :) = VS1_hrf_tmp;
end
VS_hrf = VS1_hrf(:, 1:3465);

%VS_hrf = CVS_hrf + FVS_hrf;

%% action selection(AS)
AS = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:size(action_time,1)
        AS(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1; %repelem(1,10);
    end 
end

for si = 1:Nsubj_neuro - length(rm_sub2)
    AS_hrf_tmp = conv(AS(si,:),v);
    AS_hrf(si, :) = AS_hrf_tmp;
end
AS_hrf = AS_hrf(:, 1:3465);

%% example plot
close all

sub_id = 9;
figure
subplot(5,1,1), area(1:length(Precision(sub_id,:)), Precision(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Policy Precision');t.FontSize = 8;

subplot(5,1,2), area(1:length(Pos_prediction_error(sub_id,:)), Pos_prediction_error(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Positive Predction Error');t.FontSize = 8;

subplot(5,1,3), area(1:length(Neg_prediction_error(sub_id,:)), Neg_prediction_error(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Negative Predction Error');t.FontSize = 8;

subplot(5,1,4), area(1:length(FVS_w(sub_id,:)), FVS_w(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Win outcome');t.FontSize = 8;

subplot(5,1,5), area(1:length(FVS_l(sub_id,:)), FVS_l(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Lose outcome');t.FontSize = 8;

set(gcf, 'color', 'w');
exportgraphics(gcf, [root_path, 'results/Figures/Fig_2d1.pdf'])%,'BackgroundColor','none','ContentType','vector')

figure;
subplot(5,1,1), plot(1:length(Precision_hrf(sub_id,:)), Precision_hrf(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Policy Precision hrf'); t.FontSize = 8;

subplot(5,1,2), plot(1:length(Pos_prediction_error_hrf(sub_id,:)), Pos_prediction_error_hrf(sub_id,:))
ax = gca; ax.FontSize = 6;
t = title('Positive prediction error hrf');t.FontSize = 8; 

subplot(5,1,3), plot(1:length(Neg_prediction_error_hrf(sub_id,:)), Neg_prediction_error_hrf(sub_id,:)); 
ax = gca; ax.FontSize = 6;
t = title('Negative prediction error hrf');t.FontSize = 8;

subplot(5,1,4), plot(1:length(FVS_w_hrf(sub_id,:)), FVS_w_hrf(sub_id,:)); 
ax = gca; ax.FontSize = 6;
t = title('Win outcome hrf');t.FontSize = 8;

subplot(5,1,5), plot(1:length(FVS_l_hrf(sub_id,:)), FVS_l_hrf(sub_id,:)); 
ax = gca; ax.FontSize = 6;
t = title('Lose outcome hrf');t.FontSize = 8;

%sgtitle(['Subject = ', num2str(sub_id)])
set(gcf, 'color', 'w');
exportgraphics(gcf, [root_path, 'results/Figures/Fig_2d2.pdf'])%,'BackgroundColor','none','ContentType','vector')

%%
sub_id = 9;
f1 = figure;
subplot(5,2,1), area(1:length(ST(sub_id,:)), ST(sub_id,:))
title('ST')
subplot(5,2,3), area(1:length(AS(sub_id,:)), AS(sub_id,:))
title('AS')
subplot(5,2,5), area(1:length(CVS(sub_id,:)), CVS(sub_id,:))
title('CVS')
subplot(5,2,7), area(1:length(FVS_w(sub_id,:)), FVS_w(sub_id,:))
title('FVS w')
subplot(5,2,9), area(1:length(FVS_l(sub_id,:)), FVS_l(sub_id,:))
title('FVS l')
%subplot(6,2,11), area(1:length(VS(sub_id,:)), VS(sub_id,:))
%title('VS')

subplot(5, 2, 2), plot(x(1:size(ST_hrf,2)), ST_hrf(sub_id,:), 'k'); title('ST hrf')
subplot(5, 2, 4), plot(x(1:size(AS_hrf,2)), AS_hrf(sub_id,:), 'k'); title('AS hrf')
subplot(5, 2, 6), plot(x(1:size(CVS_hrf,2)), CVS_hrf(sub_id,:), 'k'); title('CVS hrf')
subplot(5, 2, 8), plot(x(1:size(FVS_w_hrf,2)), FVS_w_hrf(sub_id,:), 'k'); title('FVS w hrf')
subplot(5, 2, 10), plot(x(1:size(FVS_l_hrf,2)), FVS_l_hrf(sub_id,:), 'k'); title('FVS l hrf')
%subplot(6, 2, 12), plot(x(1:size(VS_hrf,2)), VS_hrf(sub_id,:), 'k'); title('VS hrf')

sgtitle(['Subject = ', num2str(sub_id)])
set(gcf, 'color', 'w');

%%
save([root_path, 'data/behav/signal_data.mat'], 'ST', 'VS', 'CVS', 'FVS', 'FVS_w', 'FVS_l', 'AS', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'ST_hrf', 'VS_hrf', 'CVS_hrf', 'FVS_hrf', 'FVS_w_hrf', 'FVS_l_hrf', 'AS_hrf', '-append')

%%
PW = Pos_prediction_error .* FVS_w;
NW = Neg_prediction_error .* FVS_w;
PL = Pos_prediction_error .* FVS_l;
NL = Neg_prediction_error .* FVS_l;
PW_hrf = [];
NW_hrf = [];
PL_hrf = [];
NL_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    PW_hrf_tmp = conv(PW(si,:),v);
    PW_hrf(si, :) = PW_hrf_tmp;
    
    NW_hrf_tmp = conv(NW(si,:),v);
    NW_hrf(si, :) = NW_hrf_tmp;   
    
    PL_hrf_tmp = conv(PL(si,:),v);
    PL_hrf(si, :) = PL_hrf_tmp;
    
    NL_hrf_tmp = conv(NL(si,:),v);
    NL_hrf(si, :) = NL_hrf_tmp;
end
PW_hrf = PW_hrf(:, 1:3465);
NW_hrf = NW_hrf(:, 1:3465);
PL_hrf = PL_hrf(:, 1:3465);
NL_hrf = NL_hrf(:, 1:3465);
%%
figure;
subplot(4,2,1), area(1:length(PW(sub_id,:)), PW(sub_id,:))
title('PW')
subplot(4,2,3), area(1:length(NW(sub_id,:)), NW(sub_id,:))
title('NW')
subplot(4,2,5), area(1:length(PL(sub_id,:)), PL(sub_id,:))
title('PL')
subplot(4,2,7), area(1:length(NL(sub_id,:)), NL(sub_id,:))
title('NL')

subplot(4, 2, 2), plot(x(1:size(PW_hrf,2)), PW_hrf(sub_id,:)); title('PW hrf')
subplot(4, 2, 4), plot(x(1:size(NW_hrf,2)), NW_hrf(sub_id,:)); title('NW hrf')
subplot(4, 2, 6), plot(x(1:size(PL_hrf,2)), PL_hrf(sub_id,:)); title('PL hrf')
subplot(4, 2, 8), plot(x(1:size(NL_hrf,2)), NL_hrf(sub_id,:)); title('NL hrf')
set(gcf, 'color', 'w');
save([root_path, 'data/behav/signal_data.mat'], 'PW', 'NW', 'PL','NL', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'PW_hrf', 'NW_hrf', 'PL_hrf', 'NL_hrf', '-append')

%%
load([root_path,'data/behav/feedtype_ts.mat'], 'loss_ts', 'reward_ts')
loss_incong_ts = loss_ts(:,2) + loss_ts(:,3);
reward_incong_ts = reward_ts(:,2) + reward_ts(:,3);
uncertain_trials = [find(loss_incong_ts > 0.2); find(reward_incong_ts > 0.2) + 40];
uncertain_trials_1 = uncertain_trials - 1;
uncertain_trials = uncertain_trials_1;
uncertain_trials(find(uncertain_trials <= 0)) = [];
uncertain_trials(find(uncertain_trials == 41-1)) = [];
save([root_path,'data/behav/feedtype_ts.mat'], 'uncertain_trials', '-append')

uncertain_trials_1 = []; uncertain_trials_2 = [];
for ii = 1:length(uncertain_trials)
    if uncertain_trials(ii) < 41
        if uncertain_trials(ii) < 10
            uncertain_trials_1 = [uncertain_trials_1; uncertain_trials(ii)];
        else
            uncertain_trials_2 = [uncertain_trials_2; uncertain_trials(ii)];
        end
    else
        if uncertain_trials(ii) < 50
            uncertain_trials_1 = [uncertain_trials_1; uncertain_trials(ii)];
        else
            uncertain_trials_2 = [uncertain_trials_2; uncertain_trials(ii)];
        end
    end
end
save([root_path,'data/behav/feedtype_ts.mat'], 'uncertain_trials_1', 'uncertain_trials_2', '-append')

FVS_uncertain = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_certain = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(feed_time)
        if ismember(ti, uncertain_trials)
            FVS_uncertain(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        else
            FVS_certain(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        end
    end
end
FVS_uncertain_hrf = []; FVS_certain_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    FVS_uncertain_hrf(si, :) = conv(FVS_uncertain(si,:),v);
    FVS_certain_hrf(si, :) = conv(FVS_certain(si,:),v);
end
FVS_uncertain_hrf = FVS_uncertain_hrf(:, 1:3465);
FVS_certain_hrf = FVS_certain_hrf(:, 1:3465);
%
FVS_uncertain1 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_uncertain2 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(feed_time)
        if ismember(ti, uncertain_trials_1)
            FVS_uncertain1(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2)
            FVS_uncertain2(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1;
        end
    end
end
FVS_uncertain1_hrf = []; FVS_uncertain2_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    FVS_uncertain1_hrf(si, :) = conv(FVS_uncertain1(si,:),v);
    FVS_uncertain2_hrf(si, :) = conv(FVS_uncertain2(si,:),v);
end
FVS_uncertain1_hrf = FVS_uncertain1_hrf(:, 1:3465);
FVS_uncertain2_hrf = FVS_uncertain2_hrf(:, 1:3465);

%
PE_uncertain = UnSigned_RPE .* FVS_uncertain;
PE_certain = UnSigned_RPE .* FVS_certain;
PE_uncertain_hrf = [];
PE_certain_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    PE_uncertain_hrf(si, :) = conv(PE_uncertain(si,:),v);
    PE_certain_hrf(si, :) = conv(PE_certain(si,:),v);
end
PE_uncertain_hrf = PE_uncertain_hrf(:, 1:3465);
PE_certain_hrf = PE_certain_hrf(:, 1:3465);

figure;
subplot(6, 2, 1), area(1:length(PE_certain(sub_id,:)), PE_certain(sub_id,:))
title('PE certain')
subplot(6, 2, 3), area(1:length(PE_uncertain(sub_id,:)), PE_uncertain(sub_id,:))
title('PE uncertain')
subplot(6, 2, 5), area(1:length(FVS_certain(sub_id,:)), FVS_certain(sub_id,:))
title('FVS certain')
subplot(6, 2, 7), area(1:length(FVS_uncertain(sub_id,:)), FVS_uncertain(sub_id,:))
title('FVS uncertain')
subplot(6, 2, 9), area(1:length(FVS_uncertain1(sub_id,:)), FVS_uncertain1(sub_id,:))
title('FVS uncertain1')
subplot(6, 2, 11), area(1:length(FVS_uncertain2(sub_id,:)), FVS_uncertain2(sub_id,:))
title('FVS uncertain2')

subplot(6, 2, 2), plot(x(1:size(PE_certain_hrf,2)), PE_certain_hrf(sub_id,:)); title('PE certain hrf')
subplot(6, 2, 4), plot(x(1:size(PE_uncertain_hrf,2)), PE_uncertain_hrf(sub_id,:)); title('PE uncertain hrf')
subplot(6, 2, 6), plot(x(1:size(FVS_certain_hrf,2)), FVS_certain_hrf(sub_id,:)); title('FVS certain hrf')
subplot(6, 2, 8), plot(x(1:size(FVS_uncertain_hrf,2)), FVS_uncertain_hrf(sub_id,:)); title('FVS uncertain hrf')
subplot(6, 2, 10), plot(x(1:size(FVS_uncertain1_hrf,2)), FVS_uncertain1_hrf(sub_id,:)); title('FVS uncertain1 hrf')
subplot(6, 2, 12), plot(x(1:size(FVS_uncertain2_hrf,2)), FVS_uncertain2_hrf(sub_id,:)); title('FVS uncertain2 hrf')
set(gcf, 'color', 'w');

save([root_path, 'data/behav/signal_data.mat'], 'FVS_uncertain', 'FVS_uncertain1', 'FVS_uncertain2', 'FVS_certain', 'PE_certain', 'PE_uncertain', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'FVS_uncertain_hrf', 'FVS_uncertain1_hrf', 'FVS_uncertain2_hrf', 'FVS_certain_hrf', 'PE_certain_hrf', 'PE_uncertain_hrf', '-append')
%%
uncertain_trials_2_1 = 15; uncertain_trials_2_2 = 23;
uncertain_trials_2_3 = 62; uncertain_trials_2_4 = 76; uncertain_trials_2_5 = 78;
save([root_path,'data/behav/feedtype_ts.mat'], 'uncertain_trials_2_1', 'uncertain_trials_2_2', 'uncertain_trials_2_3', 'uncertain_trials_2_4', 'uncertain_trials_2_5', '-append')

FVS_uncertain2_1 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_uncertain2_2 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_uncertain2_3 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_uncertain2_4 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
FVS_uncertain2_5 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(feed_time)
        if ismember(ti, uncertain_trials_2_1)
            FVS_uncertain2_1(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_2)
            FVS_uncertain2_2(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1;
        elseif ismember(ti, uncertain_trials_2_3)
            FVS_uncertain2_3(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_4)
            FVS_uncertain2_4(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_5)
            FVS_uncertain2_5(si, round(feed_time(ti)*10):round(feed_time(ti)*10) + 9) = 1; 
        end
    end
end
FVS_uncertain2_1_hrf = []; FVS_uncertain2_2_hrf = []; FVS_uncertain2_3_hrf = []; FVS_uncertain2_4_hrf = []; FVS_uncertain2_5_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    FVS_uncertain2_1_hrf(si, :) = conv(FVS_uncertain2_1(si,:),v);
    FVS_uncertain2_2_hrf(si, :) = conv(FVS_uncertain2_2(si,:),v);
    FVS_uncertain2_3_hrf(si, :) = conv(FVS_uncertain2_3(si,:),v);
    FVS_uncertain2_4_hrf(si, :) = conv(FVS_uncertain2_4(si,:),v);
    FVS_uncertain2_5_hrf(si, :) = conv(FVS_uncertain2_5(si,:),v);
end
FVS_uncertain2_1_hrf = FVS_uncertain2_1_hrf(:, 1:3465);
FVS_uncertain2_2_hrf = FVS_uncertain2_2_hrf(:, 1:3465);
FVS_uncertain2_3_hrf = FVS_uncertain2_3_hrf(:, 1:3465);
FVS_uncertain2_4_hrf = FVS_uncertain2_4_hrf(:, 1:3465);
FVS_uncertain2_5_hrf = FVS_uncertain2_5_hrf(:, 1:3465);

figure;
subplot(5, 2, 1), area(1:length(FVS_uncertain2_1(sub_id,:)), FVS_uncertain2_1(sub_id,:))
subplot(5, 2, 3), area(1:length(FVS_uncertain2_2(sub_id,:)), FVS_uncertain2_2(sub_id,:))
subplot(5, 2, 5), area(1:length(FVS_uncertain2_3(sub_id,:)), FVS_uncertain2_3(sub_id,:))
subplot(5, 2, 7), area(1:length(FVS_uncertain2_4(sub_id,:)), FVS_uncertain2_4(sub_id,:))
subplot(5, 2, 9), area(1:length(FVS_uncertain2_5(sub_id,:)), FVS_uncertain2_5(sub_id,:))

subplot(5, 2, 2), plot(x(1:size(FVS_uncertain2_1_hrf,2)), FVS_uncertain2_1_hrf(sub_id,:)); 
subplot(5, 2, 4), plot(x(1:size(FVS_uncertain2_2_hrf,2)), FVS_uncertain2_2_hrf(sub_id,:));
subplot(5, 2, 6), plot(x(1:size(FVS_uncertain2_3_hrf,2)), FVS_uncertain2_3_hrf(sub_id,:)); 
subplot(5, 2, 8), plot(x(1:size(FVS_uncertain2_4_hrf,2)), FVS_uncertain2_4_hrf(sub_id,:));
subplot(5, 2, 10), plot(x(1:size(FVS_uncertain2_5_hrf,2)), FVS_uncertain2_5_hrf(sub_id,:));
set(gcf, 'color', 'w');

save([root_path, 'data/behav/signal_data.mat'], 'FVS_uncertain2_1', 'FVS_uncertain2_2', 'FVS_uncertain2_3', 'FVS_uncertain2_4', 'FVS_uncertain2_5', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'FVS_uncertain2_1_hrf', 'FVS_uncertain2_2_hrf', 'FVS_uncertain2_3_hrf', 'FVS_uncertain2_4_hrf', 'FVS_uncertain2_5_hrf', '-append')
%%
CVS_uncertain2_1 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
CVS_uncertain2_2 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
CVS_uncertain2_3 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
CVS_uncertain2_4 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
CVS_uncertain2_5 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(cue_time)
        if ismember(ti, uncertain_trials_2_1 + 1)
            CVS_uncertain2_1(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_2 + 1)
            CVS_uncertain2_2(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1;
        elseif ismember(ti, uncertain_trials_2_3 + 1)
            CVS_uncertain2_3(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_4 + 1)
            CVS_uncertain2_4(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_5 + 1)
            CVS_uncertain2_5(si, round(cue_time(ti)*10):round(cue_time(ti)*10) + 9) = 1; 
        end
    end
end
CVS_uncertain2_1_hrf = []; CVS_uncertain2_2_hrf = []; CVS_uncertain2_3_hrf = []; CVS_uncertain2_4_hrf = []; CVS_uncertain2_5_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    CVS_uncertain2_1_hrf(si, :) = conv(CVS_uncertain2_1(si,:),v);
    CVS_uncertain2_2_hrf(si, :) = conv(CVS_uncertain2_2(si,:),v);
    CVS_uncertain2_3_hrf(si, :) = conv(CVS_uncertain2_3(si,:),v);
    CVS_uncertain2_4_hrf(si, :) = conv(CVS_uncertain2_4(si,:),v);
    CVS_uncertain2_5_hrf(si, :) = conv(CVS_uncertain2_5(si,:),v);
end
CVS_uncertain2_1_hrf = CVS_uncertain2_1_hrf(:, 1:3465);
CVS_uncertain2_2_hrf = CVS_uncertain2_2_hrf(:, 1:3465);
CVS_uncertain2_3_hrf = CVS_uncertain2_3_hrf(:, 1:3465);
CVS_uncertain2_4_hrf = CVS_uncertain2_4_hrf(:, 1:3465);
CVS_uncertain2_5_hrf = CVS_uncertain2_5_hrf(:, 1:3465);

figure;
subplot(5, 2, 1), area(1:length(CVS_uncertain2_1(sub_id,:)), CVS_uncertain2_1(sub_id,:))
subplot(5, 2, 3), area(1:length(CVS_uncertain2_2(sub_id,:)), CVS_uncertain2_2(sub_id,:))
subplot(5, 2, 5), area(1:length(CVS_uncertain2_3(sub_id,:)), CVS_uncertain2_3(sub_id,:))
subplot(5, 2, 7), area(1:length(CVS_uncertain2_4(sub_id,:)), CVS_uncertain2_4(sub_id,:))
subplot(5, 2, 9), area(1:length(CVS_uncertain2_5(sub_id,:)), CVS_uncertain2_5(sub_id,:))

subplot(5, 2, 2), plot(x(1:size(CVS_uncertain2_1_hrf,2)), CVS_uncertain2_1_hrf(sub_id,:)); 
subplot(5, 2, 4), plot(x(1:size(CVS_uncertain2_2_hrf,2)), CVS_uncertain2_2_hrf(sub_id,:));
subplot(5, 2, 6), plot(x(1:size(CVS_uncertain2_3_hrf,2)), CVS_uncertain2_3_hrf(sub_id,:)); 
subplot(5, 2, 8), plot(x(1:size(CVS_uncertain2_4_hrf,2)), CVS_uncertain2_4_hrf(sub_id,:));
subplot(5, 2, 10), plot(x(1:size(CVS_uncertain2_5_hrf,2)), CVS_uncertain2_5_hrf(sub_id,:));
set(gcf, 'color', 'w');

save([root_path, 'data/behav/signal_data.mat'], 'CVS_uncertain2_1', 'CVS_uncertain2_2', 'CVS_uncertain2_3', 'CVS_uncertain2_4', 'CVS_uncertain2_5', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'CVS_uncertain2_1_hrf', 'CVS_uncertain2_2_hrf', 'CVS_uncertain2_3_hrf', 'CVS_uncertain2_4_hrf', 'CVS_uncertain2_5_hrf', '-append')

%%
AS_uncertain2_1 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
AS_uncertain2_2 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
AS_uncertain2_3 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
AS_uncertain2_4 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
AS_uncertain2_5 = zeros(Nsubj_neuro - length(rm_sub2), tp*TR*10);
for si = 1:Nsubj_neuro - length(rm_sub2)
    for ti = 1:length(cue_time)
        if ismember(ti, uncertain_trials_2_1 + 1)
            AS_uncertain2_1(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_2 + 1)
            AS_uncertain2_2(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1;
        elseif ismember(ti, uncertain_trials_2_3 + 1)
            AS_uncertain2_3(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_4 + 1)
            AS_uncertain2_4(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1; 
        elseif ismember(ti, uncertain_trials_2_5 + 1)
            AS_uncertain2_5(si, round((cue_time(ti)+action_time(ti,si))*10):round((cue_time(ti)+action_time(ti,si))*10) + 9) = 1; 
        end
    end
end
AS_uncertain2_1_hrf = []; AS_uncertain2_2_hrf = []; AS_uncertain2_3_hrf = []; AS_uncertain2_4_hrf = []; AS_uncertain2_5_hrf = [];
for si = 1:Nsubj_neuro - length(rm_sub2)
    AS_uncertain2_1_hrf(si, :) = conv(AS_uncertain2_1(si,:),v);
    AS_uncertain2_2_hrf(si, :) = conv(AS_uncertain2_2(si,:),v);
    AS_uncertain2_3_hrf(si, :) = conv(AS_uncertain2_3(si,:),v);
    AS_uncertain2_4_hrf(si, :) = conv(AS_uncertain2_4(si,:),v);
    AS_uncertain2_5_hrf(si, :) = conv(AS_uncertain2_5(si,:),v);
end
AS_uncertain2_1_hrf = AS_uncertain2_1_hrf(:, 1:3465);
AS_uncertain2_2_hrf = AS_uncertain2_2_hrf(:, 1:3465);
AS_uncertain2_3_hrf = AS_uncertain2_3_hrf(:, 1:3465);
AS_uncertain2_4_hrf = AS_uncertain2_4_hrf(:, 1:3465);
AS_uncertain2_5_hrf = AS_uncertain2_5_hrf(:, 1:3465);

figure;
subplot(5, 2, 1), area(1:length(AS_uncertain2_1(sub_id,:)), AS_uncertain2_1(sub_id,:))
subplot(5, 2, 3), area(1:length(AS_uncertain2_2(sub_id,:)), AS_uncertain2_2(sub_id,:))
subplot(5, 2, 5), area(1:length(AS_uncertain2_3(sub_id,:)), AS_uncertain2_3(sub_id,:))
subplot(5, 2, 7), area(1:length(AS_uncertain2_4(sub_id,:)), AS_uncertain2_4(sub_id,:))
subplot(5, 2, 9), area(1:length(AS_uncertain2_5(sub_id,:)), AS_uncertain2_5(sub_id,:))

subplot(5, 2, 2), plot(x(1:size(AS_uncertain2_1_hrf,2)), AS_uncertain2_1_hrf(sub_id,:)); 
subplot(5, 2, 4), plot(x(1:size(AS_uncertain2_2_hrf,2)), AS_uncertain2_2_hrf(sub_id,:));
subplot(5, 2, 6), plot(x(1:size(AS_uncertain2_3_hrf,2)), AS_uncertain2_3_hrf(sub_id,:)); 
subplot(5, 2, 8), plot(x(1:size(AS_uncertain2_4_hrf,2)), AS_uncertain2_4_hrf(sub_id,:));
subplot(5, 2, 10), plot(x(1:size(AS_uncertain2_5_hrf,2)), AS_uncertain2_5_hrf(sub_id,:));
set(gcf, 'color', 'w');

save([root_path, 'data/behav/signal_data.mat'], 'AS_uncertain2_1', 'AS_uncertain2_2', 'AS_uncertain2_3', 'AS_uncertain2_4', 'AS_uncertain2_5', '-append')
save([root_path, 'data/behav/signal_data.mat'], 'AS_uncertain2_1_hrf', 'AS_uncertain2_2_hrf', 'AS_uncertain2_3_hrf', 'AS_uncertain2_4_hrf', 'AS_uncertain2_5_hrf', '-append')

%% 
figure;
sub_id = 15;
subplot(5, 1, 1), plot(x(1:size(FVS_uncertain2_1_hrf,2)), FVS_uncertain2_1_hrf(sub_id,:), 'color', [0.3010 0.7450 0.9330]); 
hold on; plot(x(1:size(CVS_uncertain2_1_hrf,2)), CVS_uncertain2_1_hrf(sub_id,:), 'color', [0.4660 0.6740 0.1880]); 
hold on; plot(x(1:size(AS_uncertain2_1_hrf,2)), AS_uncertain2_1_hrf(sub_id,:), 'color', [0.9290 0.6940 0.1250]); 
ylabel({'Loss task ';'Trial 15-16'}, 'HorizontalAlignment', 'right', 'Rotation', 0)
%yline(0, 'color', 'w', 'LineWidth', 1, 'Alpha', 1)
subplot(5, 1, 2), plot(x(1:size(FVS_uncertain2_2_hrf,2)), FVS_uncertain2_2_hrf(sub_id,:), 'color', [0.3010 0.7450 0.9330]);
hold on; plot(x(1:size(CVS_uncertain2_2_hrf,2)), CVS_uncertain2_2_hrf(sub_id,:), 'color', [0.4660 0.6740 0.1880]);
hold on; plot(x(1:size(AS_uncertain2_2_hrf,2)), AS_uncertain2_2_hrf(sub_id,:), 'color', [0.9290 0.6940 0.1250]);
ylabel({'Loss task ';'Trial 23-24'}, 'HorizontalAlignment', 'right', 'Rotation', 0)
subplot(5, 1, 3), plot(x(1:size(FVS_uncertain2_3_hrf,2)), FVS_uncertain2_3_hrf(sub_id,:), 'color', [0.3010 0.7450 0.9330]); 
hold on; plot(x(1:size(CVS_uncertain2_3_hrf,2)), CVS_uncertain2_3_hrf(sub_id,:), 'color', [0.4660 0.6740 0.1880]); 
hold on; plot(x(1:size(AS_uncertain2_3_hrf,2)), AS_uncertain2_3_hrf(sub_id,:), 'color', [0.9290 0.6940 0.1250]); 
ylabel({'Reward task';'Trial 22-23  '}, 'HorizontalAlignment', 'right', 'Rotation', 0)
subplot(5, 1, 4), plot(x(1:size(FVS_uncertain2_4_hrf,2)), FVS_uncertain2_4_hrf(sub_id,:), 'color', [0.3010 0.7450 0.9330]);
hold on; plot(x(1:size(CVS_uncertain2_4_hrf,2)), CVS_uncertain2_4_hrf(sub_id,:), 'color', [0.4660 0.6740 0.1880]);
hold on; plot(x(1:size(AS_uncertain2_4_hrf,2)), AS_uncertain2_4_hrf(sub_id,:), 'color', [0.9290 0.6940 0.1250]);
ylabel({'Reward task';'Trial 36-37  '}, 'HorizontalAlignment', 'right', 'Rotation', 0)
subplot(5, 1, 5), plot(x(1:size(FVS_uncertain2_5_hrf,2)), FVS_uncertain2_5_hrf(sub_id,:), 'color', [0.3010 0.7450 0.9330]);
hold on; plot(x(1:size(CVS_uncertain2_5_hrf,2)), CVS_uncertain2_5_hrf(sub_id,:), 'color', [0.4660 0.6740 0.1880]);
hold on; plot(x(1:size(AS_uncertain2_5_hrf,2)), AS_uncertain2_5_hrf(sub_id,:), 'color', [0.9290 0.6940 0.1250]);
ylabel({'Reward task';'Trial 38-39  '}, 'HorizontalAlignment', 'right', 'Rotation', 0)
set(gcf, 'color', 'w');
set(gcf, 'Position', [680 255 1101 723])
exportgraphics(gcf, [root_path, 'results/Figures/Fig_6c_plot.pdf'])%,'BackgroundColor','none','ContentType','vector')

%%
figure; plot(1:100, hrf(1:100))
function y = hrf(x)
	peak_values = pdf('Gamma', x, 6, 1);
	undershoot_values = pdf('Gamma', x, 12, 1);
	values = peak_values - 0.35 * undershoot_values;
	y = (values ./ max(values)) * 0.6;
end