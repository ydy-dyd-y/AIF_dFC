clear
clc
close all

str1 = pwd; str = split(str1, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end
addpath(genpath([root_path, 'tools/spm12/']))
addpath(genpath('/APP/enhpc/NeuroScience/NeuroscienceToolbox/tapas'))

load([root_path, 'data/behav/behavior_data.mat'])
Nsubj = numel(data);

precL1 = []; precR1 = []; varL1 = []; varR1 = []; peL1 = []; peR1 = [];
precL2 = []; precR2 = []; varL2 = []; varR2 = []; peL2 = []; peR2 = [];
precL3 = []; precR3 = []; varL3 = []; varR3 = []; peL3 = []; peR3 = [];

try load([root_path, 'results/behav/hgf.mat'])
catch est_hgf_lore = cell(Nsubj, 2);
end

for i = 1:Nsubj
    if isempty(est_hgf_lore{i,1})
        est_hgf_lore{i,1} = tapas_fitModel(data{i}.response{1}, data{i}.input{1}, 'tapas_hgf_binary_config','tapas_softmax_binary_config');
        est_hgf_lore{i,2} = tapas_fitModel(data{i}.response{2}, data{i}.input{2}, 'tapas_hgf_binary_config','tapas_softmax_binary_config');
    end
    % volatility prediction error
    peL1 = [peL1, est_hgf_lore{i,1}.traj.da(:,1)];
    peR1 = [peR1, est_hgf_lore{i,2}.traj.da(:,1)]; 
    peL2 = [peL2, est_hgf_lore{i,1}.traj.da(:,2)];
    peR2 = [peR2, est_hgf_lore{i,2}.traj.da(:,2)];     
    peL3 = [peL3, est_hgf_lore{i,1}.traj.da(:,3)];
    peR3 = [peR3, est_hgf_lore{i,2}.traj.da(:,3)]; 
    
    % variance of predictions
    varL1 = [varL1, est_hgf_lore{i,1}.traj.sahat(:,1)];
    varR1 = [varR1, est_hgf_lore{i,2}.traj.sahat(:,1)];    
    varL2 = [varL2, est_hgf_lore{i,1}.traj.sahat(:,2)];
    varR2 = [varR2, est_hgf_lore{i,2}.traj.sahat(:,2)];
    varL3 = [varL3, est_hgf_lore{i,1}.traj.sahat(:,3)];
    varR3 = [varR3, est_hgf_lore{i,2}.traj.sahat(:,3)];    
    
    % precision of predictions
    precL1 = [precL1, est_hgf_lore{i,1}.traj.psi(:,1)];
    precR1 = [precR1, est_hgf_lore{i,2}.traj.psi(:,1)];    
    precL2 = [precL2, est_hgf_lore{i,1}.traj.psi(:,2)];
    precR2 = [precR2, est_hgf_lore{i,2}.traj.psi(:,2)];
    precL3 = [precL3, est_hgf_lore{i,1}.traj.psi(:,3)];
    precR3 = [precR3, est_hgf_lore{i,2}.traj.psi(:,3)];    
end

save([root_path, 'results/behav/hgf.mat'], 'est_hgf_lore')
rm_sub1 = [20, 32];
rm_sub2 = [22, 36, 44, 57];
varL1(:,rm_sub1) = []; varR1(:,rm_sub1) = [];
varL1(:,rm_sub2) = []; varR1(:,rm_sub2) = [];

%%
figure; heatmap(peL3')
title('prediction error level3')
%%
figure;
subplot(1, 2, 1);
inBetween1 = [mean(varL1, 2)' + std(varL1'), fliplr(mean(varL1, 2)' - std(varL1'))];
inBetween2 = [max(varL1'), fliplr(mean(varL1, 2)' - std(varL1'))]; %fliplr(min(varL1'))];

fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(varL1, 2), 'k', 'LineWidth', 2)
ax = gca;
ax.FontSize = 6;
%{
subplot(2, 2, 3);
inBetween1 = [mean(abs(peL1), 2)' + std(abs(peL1)'), fliplr(mean(abs(peL1), 2)' - std(abs(peL1)'))];
inBetween2 = [max(peL1'), fliplr(mean(abs(peL1), 2)' - std(abs(peL1)'))]; %fliplr(min(peL1'))];
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(abs(peL1), 2), 'r', 'LineWidth', 2)
%}
subplot(1, 2, 2); 
inBetween1 = [mean(varR1, 2)' + std(varR1'), fliplr(mean(varR1, 2)' - std(varR1'))];
inBetween2 = [max(varR1'), fliplr(mean(varR1, 2)' - std(varR1'))]; %fliplr(min(varR1'))];
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(varR1, 2), 'k', 'LineWidth', 2)
ax = gca;
ax.FontSize = 6;
%{
subplot(2, 2, 4);
inBetween1 = [mean(abs(peR1), 2)' + std(abs(peR1)'), fliplr(mean(abs(peR1), 2)' - std(abs(peR1)'))];
inBetween2 = [max(peR1'), fliplr(mean(abs(peR1), 2)' - std(abs(peR1)'))]; %fliplr(min(peR1'))];
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(abs(peR1), 2), 'r', 'LineWidth', 2)
hold on
plot(1:40, mean(abs(peR1), 2), 'r', 'LineWidth', 2)
%}
sgtitle('Inverse of prediction precision')
set(gcf, 'color', 'w')
f = gcf;
f.Position = [105 79 1581 374];
exportgraphics(gcf, [root_path, 'results/Figures/figure5_supp.pdf'])%,'BackgroundColor','none','ContentType','vector')
%%
figure;
subplot(1, 2, 1);
for si = 1:size(varL1, 2)
end
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(varL1, 2), 'k', 'LineWidth', 2)
%{
subplot(2, 2, 3);
inBetween1 = [mean(abs(peL1), 2)' + std(abs(peL1)'), fliplr(mean(abs(peL1), 2)' - std(abs(peL1)'))];
inBetween2 = [max(peL1'), fliplr(mean(abs(peL1), 2)' - std(abs(peL1)'))]; %fliplr(min(peL1'))];
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(abs(peL1), 2), 'r', 'LineWidth', 2)
%}
subplot(1, 2, 2); 
inBetween1 = [mean(varR1, 2)' + std(varR1'), fliplr(mean(varR1, 2)' - std(varR1'))];
inBetween2 = [max(varR1'), fliplr(mean(varR1, 2)' - std(varR1'))]; %fliplr(min(varR1'))];
fill([1:40, fliplr(1:40)], inBetween2, 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
hold on
plot(1:40, mean(varR1, 2), 'k', 'LineWidth', 2)
%%
corrL = []; corrR = [];
for si = 1:Nsubj
    tmp = corrcoef(abs(peR1(1:39,si)), precR1(2:end,si));
    corrR = [corrR, tmp(1,2)];
    tmp = corrcoef(abs(peL1(1:39,si)), precL1(2:end,si));
    corrL = [corrL, tmp(1,2)];
end
mean(corrR)
mean(corrL)

% only with  PPE trial
corrL = []; corrR = [];
for si = 1:Nsubj
    tmp = corrcoef(abs(peL1(find(data{1,1}.input{1,1} == 1),si)), precL1(find(data{1,1}.input{1,1} == 1)+1,si));
    corrL = [corrL, tmp(1,2)];
    R_idx = find(data{1,1}.input{1,2} == 1); R_idx(end) = [];
    tmp = corrcoef(abs(peR1(R_idx,si)), precR1(R_idx+1,si));
    corrR = [corrR, tmp(1,2)];
end
mean(corrR)
mean(corrL)

% only with  NPE trial
corrL = []; corrR = [];
for si = 1:Nsubj
    L_idx = find(data{1,1}.input{1,1} == 0); L_idx(end) = [];
    tmp = corrcoef(abs(peL1(L_idx,si)), precL1(L_idx+1,si));
    corrL = [corrL, tmp(1,2)];
    R_idx = find(data{1,1}.input{1,2} == 0); %R_idx(end) = [];
    tmp = corrcoef(abs(peR1(R_idx,si)), precR1(R_idx+1,si));
    corrR = [corrR, tmp(1,2)];
end
mean(corrR)
mean(corrL)
%%
figure;
n_term = 10;
for si = 1:10
    subplot(10, 1, si)
    plot(peR3(:,si+n_term))
end
sgtitle('prediction error of level3')

%%
load([root_path, 'data/behav/feedback_4type.mat'], 'trial_sub_4type_idx_loss', 'trial_sub_4type_idx_rew')  % (40 x 54)
mean_varL_type = [];
for si = 1:size(varL1,2)
    for type_i = 1:max(max(trial_sub_4type_idx_loss))
        mean_varL_type(type_i, si) = nanmean(varL1(find(trial_sub_4type_idx_loss(:,si) == type_i), si));
    end
end

mean_varR_type = [];
for si = 1:size(varR1,2)
    for type_i = 1:max(max(trial_sub_4type_idx_rew))
        mean_varR_type(type_i, si) = nanmean(varR1(find(trial_sub_4type_idx_rew(:,si) == type_i), si));
    end
end

% ttest
H_L = []; P_L = [];
[h, p] = ttest(mean_varL_type(1, :), mean_varL_type(2, :), 'Tail', 'left');
H_L = [H_L; h]; P_L = [P_L; p];
[h, p] = ttest(mean_varL_type(1, :), mean_varL_type(3, :), 'Tail', 'left');
H_L = [H_L; h]; P_L = [P_L; p];
[h, p] = ttest(mean_varL_type(4, :), mean_varL_type(2, :), 'Tail', 'left');
H_L = [H_L; h]; P_L = [P_L; p];
[h, p] = ttest(mean_varL_type(4, :), mean_varL_type(3, :), 'Tail', 'left');
H_L = [H_L; h]; P_L = [P_L; p];

H_R = []; P_R = [];
[h, p] = ttest(mean_varR_type(1, :), mean_varR_type(2, :), 'Tail', 'left');
H_R = [H_R; h]; P_R = [P_R; p];
[h, p] = ttest(mean_varR_type(1, :), mean_varR_type(3, :), 'Tail', 'left');
H_R = [H_R; h]; P_R = [P_R; p];
[h, p] = ttest(mean_varR_type(4, :), mean_varR_type(2, :), 'Tail', 'left');
H_R = [H_R; h]; P_R = [P_R; p];
[h, p] = ttest(mean_varR_type(4, :), mean_varR_type(3, :), 'Tail', 'left');
H_R = [H_R; h]; P_R = [P_R; p];

save([root_path, 'data/behav/stateVariance.mat'], 'mean_varL_type', 'mean_varR_type')
%%
load([root_path,'data/behav/feedtype_ts.mat'], 'incong_trials')
cong_trials_L = setdiff(1:40', incong_trials);
cong_trials_R = setdiff(41:80', incong_trials);

mean_varL_type2 = [];
for si = 1:size(varL1,2)
    mean_varL_type2(1, si) = nanmean(varL1(incong_trials(1:max(find(incong_trials < 41))), si));
    mean_varL_type2(2, si) = nanmean(varL1(cong_trials_L, si));    
end

mean_varR_type2 = [];
for si = 1:size(varR1,2)
    mean_varR_type2(1, si) = nanmean(varR1(incong_trials(min(find(incong_trials > 40)):end) - 40, si));
    mean_varR_type2(2, si) = nanmean(varR1(cong_trials_R - 40, si));    
end
save([root_path, 'data/behav/stateVariance.mat'], 'mean_varL_type2', 'mean_varR_type2', '-append')

%%
figure; subplot(3,1,1); heatmap(precL1')
subplot(3,1,2); heatmap(precL2')
subplot(3,1,3); heatmap(precL3')

%%
sub_id = 1;
figure; subplot(3, 1, 1); plot(varL3(:,sub_id))
subplot(3, 1, 2); plot(varL2(:,sub_id))
subplot(3, 1, 3); plot(varL1(:,sub_id))
%%
wrongL16 = zeros(Nsubj, 1);
wrongL24 = zeros(Nsubj, 1);
wrongR23 = zeros(Nsubj, 1);
wrongR37 = zeros(Nsubj, 1);
for si = 1:Nsubj
    if est_hgf_lore{si, 1}.y(16) == 0
        wrongL16(si) = 1;
    end
    if est_hgf_lore{si, 1}.y(24) == 0
        wrongL24(si) = 1;
    end
    if est_hgf_lore{si, 2}.y(23) == 0
        wrongR23(si) = 1;
    end
    if est_hgf_lore{si, 2}.y(37) == 0
        wrongR37(si) = 1;
    end
end
wrongL16(rm_sub1) = []; wrongL24(rm_sub1) = []; wrongR23(rm_sub1) = []; wrongR37(rm_sub1) = [];
wrongL16(rm_sub2) = []; wrongL24(rm_sub2) = []; wrongR23(rm_sub2) = []; wrongR37(rm_sub2) = [];

rightL16 = find(wrongL16 == 0);
wrongL16 = find(wrongL16);

rightL24 = find(wrongL24 == 0);
wrongL24 = find(wrongL24);

rightR23 = find(wrongR23 == 0);
wrongR23 = find(wrongR23);

rightR37 = find(wrongR37 == 0);
wrongR37 = find(wrongR37);

save([root_path, 'data/behav/PL_idx.mat'], 'wrongL16', 'wrongL24', 'wrongR23', 'wrongR37', 'rightL16', 'rightL24', 'rightR23', 'rightR37')
%%
for sub_id = wrongL16'
    tapas_hgf_binary_plotTraj(est_hgf_lore{sub_id, 1});
end

close all
for sub_id = wrongL24'
    tapas_hgf_binary_plotTraj(est_hgf_lore{sub_id, 1});
end

close all
for sub_id = wrongR23'
    tapas_hgf_binary_plotTraj(est_hgf_lore{sub_id, 2});
end

close all
for sub_id = wrongR37'
    tapas_hgf_binary_plotTraj(est_hgf_lore{sub_id, 2});
end

%%
muhatL16_wrong = [];
for sub_id = wrongL16'
    muhatL16_wrong = [muhatL16_wrong, est_hgf_lore{sub_id, 1}.traj.muhat(16,3)];
end
muhatL16_right = [];
for sub_id = rightL16'
    muhatL16_right = [muhatL16_right, est_hgf_lore{sub_id, 1}.traj.muhat(16,3)];
end

muhatL24_wrong = [];
for sub_id = wrongL24'
    muhatL24_wrong = [muhatL24_wrong, est_hgf_lore{sub_id, 1}.traj.muhat(24,3)];
end
muhatL24_right = [];
for sub_id = rightL24'
    muhatL24_right = [muhatL24_right, est_hgf_lore{sub_id, 1}.traj.muhat(24,3)];
end

muhatR23_wrong = [];
for sub_id = wrongR23'
    muhatR23_wrong = [muhatR23_wrong, est_hgf_lore{sub_id, 2}.traj.muhat(23,3)];
end
muhatR23_right = [];
for sub_id = rightR23'
    muhatR23_right = [muhatR23_right, est_hgf_lore{sub_id, 2}.traj.muhat(23,3)];
end

muhatR37_wrong = [];
for sub_id = wrongR37'
    muhatR37_wrong = [muhatR37_wrong, est_hgf_lore{sub_id, 2}.traj.muhat(37,3)];
end
muhatR37_right = [];
for sub_id = rightR37'
    muhatR37_right = [muhatR37_right, est_hgf_lore{sub_id, 2}.traj.muhat(37,3)];
end

mean(muhatL24_wrong)
mean(muhatL24_right)
std(muhatL24_wrong)
std(muhatL24_right)

mean(muhatR37_wrong)
mean(muhatR37_right)
std(muhatR37_wrong)
std(muhatR37_right)
%%
window = 5;
dISC = [];
for si = 1:size(varL1,2)
    x = varL1(:,si);
    tmp = varL1; tmp(:,si) = []; y = mean(tmp,2);  % (40 x 1)
    for wi = 1:size(varL1,1)-window
        x_tmp = x(wi:wi+window);
        y_tmp = y(wi:wi+window);
        isc = sum(x_tmp .* y_tmp)/sqrt(sum((x_tmp - mean(x_tmp)).^2) * sum((y_tmp - mean(y_tmp)).^2));
        dISC(si, wi) = isc;
    end
end
figure; plot(mean(dISC,1))
%% model estimation parameter
model_para_hgf_lore = zeros(110,6);

for i = N
    if isempty(data{i})
        continue
    end
    model_para_hgf_lore(i, 1) = est_hgf_lore{i, 1}.optim.LME;
    model_para_hgf_lore(i, 2) = est_hgf_lore{i, 1}.optim.AIC;
    model_para_hgf_lore(i, 3) = est_hgf_lore{i, 1}.optim.BIC;
    
    model_para_hgf_lore(i, 4) = est_hgf_lore{i, 2}.optim.LME;
    model_para_hgf_lore(i, 5) = est_hgf_lore{i, 2}.optim.AIC;
    model_para_hgf_lore(i, 6) = est_hgf_lore{i, 2}.optim.BIC;

end

%l_mean = nanmean(model_para_l_1030);
%r_mean = nanmean(model_para_r_1030);
para_hgf_lore_mean = sum(model_para_hgf_lore) / 60
