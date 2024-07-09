clear
clc
close all

str1 = pwd; str = split(str1, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

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

%% Precision
Nsubj_neuro = Nsubj - 2;
action_time(:,[20, 32]) = [];

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

%%
rm_sub = [22, 36, 44, 57];
V_54 = V;
V_54(:,:,rm_sub) = [];

V_pos_54 = V_54;
V_pos_54(find(V_pos_54 < 0)) = 0;
V_neg_54 = V_54;
V_neg_54(find(V_neg_54 > 0)) = 0;
V_neg_54 = abs(V_neg_54);

W_54 = W_dn;
%W_54 = W_wn;

W_54(:,:,rm_sub) = [];

V_pos = [];
V_neg = [];
W = [];
for ti = 1: 80
    for si = 1:54
        V_pos(ti, si) = max(V_pos_54(:,ti, si));
        V_neg(ti, si) = max(V_neg_54(:,ti, si));
        W(ti, si) = max(W_54(:,ti, si));
    end
end

feedback_win_or_lose_54 = feedback_win_or_lose;
feedback_win_or_lose_54(:, rm_sub) = [];  % (80 x 54)

%%
trial_sub_4type_idx_loss = [];
prec_4types_lose = cell(1,4);
for ti = 1:40
    for si = 1:54
        if feedback_win_or_lose_54(ti,si) == 1 && V_pos(ti,si) > 0  % prec high
            trial_sub_4type_idx_loss(ti, si) = 1;
            prec_4types_lose{1, 1} = [prec_4types_lose{1, 1}; W(ti,si)];
        elseif feedback_win_or_lose_54(ti,si) == 1 && V_neg(ti,si) > 0  % prec low
            trial_sub_4type_idx_loss(ti, si) = 2;
            prec_4types_lose{1, 2} = [prec_4types_lose{1, 2}; W(ti,si)];
        elseif feedback_win_or_lose_54(ti,si) == 0 && V_pos(ti,si) > 0  % prec low
            trial_sub_4type_idx_loss(ti, si) = 3;
            prec_4types_lose{1, 3} = [prec_4types_lose{1, 3}; W(ti,si)];
        elseif feedback_win_or_lose_54(ti,si) == 0 && V_neg(ti,si) > 0  % prec high
            trial_sub_4type_idx_loss(ti, si) = 4;
            prec_4types_lose{1, 4} = [prec_4types_lose{1, 4}; W(ti,si)];
        end
    end
end

trial_sub_4type_idx_rew = [];
prec_4types_rew = cell(1,4);
for ti = 1:40
    for si = 1:54
        if feedback_win_or_lose_54(40+ti,si) == 1 && V_pos(40+ti,si) > 0  % prec high
            trial_sub_4type_idx_rew(ti, si) = 1;
            prec_4types_rew{1, 1} = [prec_4types_rew{1, 1}; W(40+ti,si)];
        elseif feedback_win_or_lose_54(40+ti,si) == 1 && V_neg(40+ti,si) > 0  % prec low
            trial_sub_4type_idx_rew(ti, si) = 2;
            prec_4types_rew{1, 2} = [prec_4types_rew{1, 2}; W(40+ti,si)];
        elseif feedback_win_or_lose_54(40+ti,si) == 0 && V_pos(40+ti,si) > 0  % prec low
            trial_sub_4type_idx_rew(ti, si) = 3;
            prec_4types_rew{1, 3} = [prec_4types_rew{1, 3}; W(40+ti,si)];
        elseif feedback_win_or_lose_54(40+ti,si) == 0 && V_neg(40+ti,si) > 0  % prec high
            trial_sub_4type_idx_rew(ti, si) = 4;
            prec_4types_rew{1, 4} = [prec_4types_rew{1, 4}; W(40+ti,si)];
        end
    end
end

mean_prec_4types_lose = cell(1,4);
mean_prec_4types_rew = cell(1,4);

prec_4types_trialNum = cell(2,4);
for si = 1:54
    subj_wp_lose = [];
    subj_wn_lose = [];
    subj_lp_lose = [];
    subj_ln_lose = [];
    for ti = 1:39
        if feedback_win_or_lose_54(ti,si) == 1 && V_pos(ti,si) > 0  % prec high
            subj_wp_lose = [subj_wp_lose, W(ti+1,si)];
        elseif feedback_win_or_lose_54(ti,si) == 1 && V_neg(ti,si) > 0  % prec low
            subj_wn_lose = [subj_wn_lose, W(ti+1,si)];
        elseif feedback_win_or_lose_54(ti,si) == 0 && V_pos(ti,si) > 0  % prec low
            subj_lp_lose = [subj_lp_lose, W(ti+1,si)];
        elseif feedback_win_or_lose_54(ti,si) == 0 && V_neg(ti,si) > 0  % prec high
            subj_ln_lose = [subj_ln_lose, W(ti+1,si)];
        end
    end
    
    subj_wp_rew = [];
    subj_wn_rew = [];
    subj_lp_rew = [];
    subj_ln_rew = [];
    for ti = 1:39
        if feedback_win_or_lose_54(ti+40,si) == 1 && V_pos(ti+40,si) > 0  % prec high
            subj_wp_rew = [subj_wp_rew, W(ti+40+1,si)];
        elseif feedback_win_or_lose_54(ti+40,si) == 1 && V_neg(ti+40,si) > 0  % prec low
            subj_wn_rew = [subj_wn_rew, W(ti+40+1,si)];
        elseif feedback_win_or_lose_54(ti+40,si) == 0 && V_pos(ti+40,si) > 0  % prec low
            subj_lp_rew = [subj_lp_rew, W(ti+40+1,si)];
        elseif feedback_win_or_lose_54(ti+40,si) == 0 && V_neg(ti+40,si) > 0  % prec high
            subj_ln_rew = [subj_ln_rew, W(ti+40+1,si)];
        end
    end
    mean_prec_4types_lose{1, 1} = [mean_prec_4types_lose{1, 1}; mean(subj_wp_lose)];
    mean_prec_4types_lose{1, 2} = [mean_prec_4types_lose{1, 2}; mean(subj_wn_lose)];
    mean_prec_4types_lose{1, 3} = [mean_prec_4types_lose{1, 3}; mean(subj_lp_lose)];
    mean_prec_4types_lose{1, 4} = [mean_prec_4types_lose{1, 4}; mean(subj_ln_lose)];    
    prec_4types_trialNum{1, 1} = [prec_4types_trialNum{1, 1}; length(subj_wp_lose)];
    prec_4types_trialNum{1, 2} = [prec_4types_trialNum{1, 2}; length(subj_wn_lose)];
    prec_4types_trialNum{1, 3} = [prec_4types_trialNum{1, 3}; length(subj_lp_lose)];
    prec_4types_trialNum{1, 4} = [prec_4types_trialNum{1, 4}; length(subj_ln_lose)];   

    mean_prec_4types_rew{1, 1} = [mean_prec_4types_rew{1, 1}; mean(subj_wp_rew)];
    mean_prec_4types_rew{1, 2} = [mean_prec_4types_rew{1, 2}; mean(subj_wn_rew)];
    mean_prec_4types_rew{1, 3} = [mean_prec_4types_rew{1, 3}; mean(subj_lp_rew)];
    mean_prec_4types_rew{1, 4} = [mean_prec_4types_rew{1, 4}; mean(subj_ln_rew)];    
    prec_4types_trialNum{2, 1} = [prec_4types_trialNum{2, 1}; length(subj_wp_rew)];
    prec_4types_trialNum{2, 2} = [prec_4types_trialNum{2, 2}; length(subj_wn_rew)];
    prec_4types_trialNum{2, 3} = [prec_4types_trialNum{2, 3}; length(subj_lp_rew)];
    prec_4types_trialNum{2, 4} = [prec_4types_trialNum{2, 4}; length(subj_ln_rew)];  
end

save([root_path, 'data/behav/feedback_4type.mat'], 'prec_4types_lose', 'prec_4types_rew', 'mean_prec_4types_lose', 'mean_prec_4types_rew', 'trial_sub_4type_idx_loss', 'trial_sub_4type_idx_rew')

%%
load([root_path,'data/behav/feedtype_ts.mat'], 'incong_trials')

prec_2types_lose = cell(1,4);
for ti = 1:40
    for si = 1:54
        if ismember(ti, incong_trials)   % prec low
            prec_2types_lose{1, 1} = [prec_2types_lose{1, 1}; W(ti,si)];
        else                             % prec high
            prec_2types_lose{1, 2} = [prec_2types_lose{1, 2}; W(ti,si)];
        end
    end
end

prec_2types_rew = cell(1,4);
for ti = 1:40
    for si = 1:54
        if ismember(ti+40, incong_trials)    % prec low
            prec_2types_rew{1, 1} = [prec_2types_rew{1, 1}; W(40+ti,si)];
        else                                 % prec high
            prec_2types_rew{1, 2} = [prec_2types_rew{1, 2}; W(40+ti,si)];
        end
    end
end


mean_prec_2types_lose = cell(1,2);
mean_prec_2types_rew = cell(1,2);

for si = 1:54
    subj_incong_lose = [];
    subj_cong_lose = [];
    for ti = 1:39
        if ismember(ti, incong_trials)  % prec low
            subj_incong_lose = [subj_incong_lose, W(ti+1,si)];
        else 
            subj_cong_lose = [subj_cong_lose, W(ti+1,si)];
        end
    end
    
    subj_incong_rew = [];
    subj_cong_rew = [];
    for ti = 1:39
        if ismember(ti+40, incong_trials)  % prec low
            subj_incong_rew = [subj_incong_rew, W(ti+40+1,si)];
        else 
            subj_cong_rew = [subj_cong_rew, W(ti+40+1,si)];
        end
    end
    mean_prec_2types_lose{1, 1} = [mean_prec_2types_lose{1, 1}; mean(subj_incong_lose)];
    mean_prec_2types_lose{1, 2} = [mean_prec_2types_lose{1, 2}; mean(subj_cong_lose)];

    mean_prec_2types_rew{1, 1} = [mean_prec_2types_rew{1, 1}; mean(subj_incong_rew)];
    mean_prec_2types_rew{1, 2} = [mean_prec_2types_rew{1, 2}; mean(subj_cong_rew)];
end

save([root_path, 'data/behav/feedback_2type.mat'], 'prec_2types_lose', 'prec_2types_rew', 'mean_prec_2types_lose', 'mean_prec_2types_rew')

%%
loss_ts = [];
reward_ts = [];
for ti = 1:40
    for type_i = 1:4
        loss_ts(ti, type_i) = length(find(trial_sub_4type_idx_loss(ti,:) == type_i))./(Nsubj_neuro-length(rm_sub));
        reward_ts(ti, type_i) = length(find(trial_sub_4type_idx_rew(ti,:) == type_i))./(Nsubj_neuro-length(rm_sub));
    end
end

save([root_path, 'data/behav/feedtype_ts.mat'], 'loss_ts', 'reward_ts')

%%
count_trial_sub_4type_loss10 = [];
count_trial_sub_4type_rew10 = [];

for type_i = 1:4
    for si = 1:Nsubj_neuro - length(rm_sub)
        count_trial_sub_4type_loss10(type_i, si) = length(find(trial_sub_4type_idx_loss(1:10, si) == type_i));
        count_trial_sub_4type_rew10(type_i, si) = length(find(trial_sub_4type_idx_rew(1:10, si) == type_i));
        
        count_trial_sub_4type_loss30(type_i, si) = length(find(trial_sub_4type_idx_loss(11:end, si) == type_i));
        count_trial_sub_4type_rew30(type_i, si) = length(find(trial_sub_4type_idx_rew(11:end, si) == type_i));
    end
end
mean(count_trial_sub_4type_loss30, 2)
mean(count_trial_sub_4type_rew30, 2)
mean(count_trial_sub_4type_loss10, 2)
mean(count_trial_sub_4type_rew10, 2)
%%
function y = hrf(x)
	peak_values = pdf('Gamma', x, 6, 1);
	undershoot_values = pdf('Gamma', x, 12, 1);
	values = peak_values - 0.35 * undershoot_values;
	y = (values ./ max(values)) * 0.6;
end