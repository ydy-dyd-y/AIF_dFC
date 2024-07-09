%% Suppl. Fig. 3 : dynamic connectivity metrics calculated from occurrence probability estimated from GLMM 
clear
clc
close all
str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end
addpath(genpath([root_path, 'tools/']))

K = 4;
delta = 0.3;
iter_tmp = 5;
load([root_path, 'results/glmm/glmm_k',num2str(K),'_',num2str(delta),'.mat']);
W = W_iter{iter_tmp};
mus = squeeze(mus_iter(:,:,iter_tmp));
gamma_hats = squeeze(gamma_hats_iter(:,:,iter_tmp));
gamma_hats_tmp = reshape(gamma_hats, [tp, Nsubj, K]);

TR = 1.5;
 
iterations = size(gamma_hats_iter, 3);
gamma_hats_iter = reshape(gamma_hats_iter, [tp, Nsubj, K, iterations]);

partitions = NaN(tp, Nsubj, iterations);
for iter = 1:iterations
    for nn = 1:Nsubj
        tmp_gam = squeeze(gamma_hats_iter(:, nn, :, iter));  % [tp x K]
        [~, maxInd] = max(tmp_gam, [], 2);
        partitions(:, nn, iter) = maxInd;        
    end
end

partitions_tmp = partitions(:,:,iter_tmp)'; 

%%
% Dwell time, Fractional occupancy
P1=zeros(Nsubj,K);   % loss task
LT1=zeros(Nsubj,K);
P2=zeros(Nsubj,K);   % reward task
LT2=zeros(Nsubj,K);

for si=1:Nsubj
    % Select the time points representing this subject
    Ctime1=partitions_tmp(si, 8:114);
    Ctime2=partitions_tmp(si, 121:227);
    
    for c=1:K
        % Probability
        P1(si,c)=mean(Ctime1==c) + eps; % normalised by T
        P2(si,c)=mean(Ctime2==c) + eps;
        % Mean Lifetime
        Ctime_bin1=Ctime1==c;
        Ctime_bin2=Ctime2==c;
        
        % Detect switches in and out of this state
        a=find(diff(Ctime_bin1)==1);
        b=find(diff(Ctime_bin1)==-1);
                
        % We discard the cases where state sarts or ends ON
        if length(b)>length(a)
           b(1)=[];
        elseif length(a)>length(b)
           a(end)=[];
        elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
           b(1)=[];
           a(end)=[];
        end
        if ~isempty(a) && ~isempty(b)
            C_Durations1=b-a;
        else
            C_Durations1=0;
        end
        LT1(si,c)=mean(C_Durations1)*TR;
        %%%
        a=find(diff(Ctime_bin2)==1);
        b=find(diff(Ctime_bin2)==-1);
                
        if length(b)>length(a)
           b(1)=[];
        elseif length(a)>length(b)
           a(end)=[];
        elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
           b(1)=[];
           a(end)=[];
        end
        if ~isempty(a) && ~isempty(b)
            C_Durations2=b-a;
        else
            C_Durations2=0;
        end
        LT2(si,c)=mean(C_Durations2)*TR; 
    end
end
% plotting LT, P
pval_LT = nan(K,1);

for ik = 1 : K
    a = LT1(:,ik)';
    b = LT2(:,ik)';

    z = a-b;
    z_mean = mean(z);
    z_std = std(z);
    t_stats = sqrt(Nsubj) * z_mean / z_std;
    pval_LT(ik) = 2 * (1-tcdf(abs(t_stats), Nsubj-1));
end
diff_LT_state= find(pval_LT < 0.01);

%
pval_P = nan(K,1);

for ik = 1 : K
    a = P1(:,ik)';
    b = P2(:,ik)';

    z = a-b;
    z_mean = mean(z);
    z_std = std(z);
    t_stats = sqrt(Nsubj) * z_mean / z_std;
    pval_P(ik) = 2 * (1-tcdf(abs(t_stats), Nsubj-1));
end
diff_P_state= find(pval_P < 0.01);

x = 1:K;
mean_LT1 = mean(LT1)';
mean_LT2 = mean(LT2)';
vals = [mean(LT1)' mean(LT2)'];

err = [std(LT1);std(LT2)];
err = err';

f = figure;
subplot(1,2,1)
b_D = bar(vals,0.3);
hold on;

for ti = 1:size(vals, 2)
    % get x positions per group
    xpos = b_D(ti).XData + b_D(ti).XOffset;
    % draw errorbar
    errorbar(xpos, vals(:,ti), err(:,ti), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);    
end

xlabel('State #')
ylabel('Duration (seconds)')
xlim auto
legend('Loss task','Reward task', 'Location','southeast')
title('Dwell time')
    
xtips1 = b_D(1).XEndPoints;
ytips1 = b_D(1).YEndPoints;
xtips2 = b_D(2).XEndPoints;
ytips2 = b_D(2).YEndPoints;

for itx = 1:length(diff_LT_state)
    
    tmp_idx = diff_LT_state(itx);
    
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.044, '*','Color',[0 0 0],'FontSize',25, 'FontWeight','bold', 'HorizontalAlignment', 'center' )
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.027, '-','Color',[0 0 0],'FontSize',60, 'HorizontalAlignment', 'center')
    p_val_str = sprintf('p-val = %0.4f', round(pval_LT(tmp_idx),4));
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.052, p_val_str,'Color',[0 0 0],'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment', 'center' )
end

mean_P1 = mean(P1)';
mean_P2 = mean(P2)';
vals = [mean_P1 mean_P2];
err = [std(P1);std(P2)];
err = err';

subplot(1,2,2)
b_F = bar(vals,0.3);
hold on;

for ti = 1:size(vals, 2)
    % get x positions per group
    xpos = b_F(ti).XData + b_F(ti).XOffset;
    % draw errorbar
    errorbar(xpos, vals(:,ti), err(:,ti), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);    
end

xlabel('State #')
ylabel('Probability')
xlim auto
legend('Loss task','Reward task', 'Location','southeast')
title('Fractional Occupancy')

xtips1 = b_F(1).XEndPoints;
ytips1 = b_F(1).YEndPoints;
xtips2 = b_F(2).XEndPoints;
ytips2 = b_F(2).YEndPoints;
    
for itx = 1:length(diff_P_state)
    
    tmp_idx = diff_P_state(itx);
    
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.040, '*','Color',[0 0 0],'FontSize',25, 'FontWeight','bold', 'HorizontalAlignment', 'center' )
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.027, '-','Color',[0 0 0],'FontSize',60, 'HorizontalAlignment', 'center')
    p_val_str = sprintf('p-val = %0.4f', round(pval_P(tmp_idx),4));
    text((xtips1(tmp_idx)+xtips2(tmp_idx))/2, max([ytips1(tmp_idx) ytips2(tmp_idx)])+0.052, p_val_str,'Color',[0 0 0],'FontSize',8, 'FontWeight','bold', 'HorizontalAlignment', 'center' )
end
set(gcf,'color','w');
set(gcf,'position',[2141 252 1423 673]);
exportgraphics(f, [root_path, 'results/Figures/Sup_fig2.pdf'])
%% Transition Probability
task_1 = [];
task_2 = [];

for si =1:Nsubj
    task_1 = [task_1; partitions_tmp(si, 7:113)];
    task_2 = [task_2; partitions_tmp(si, 120:226)];
end

subjInd = [];
for si = 1:Nsubj
    subjInd = [subjInd; ones(size(task_1, 2), 1)*si]; % [tp*Nsubj x 1]
end

task_1 = reshape(task_1, size(subjInd));
task_2 = reshape(task_2, size(subjInd));

persist = 1; % 1: persist, 0: noPersist
start = 3; % 1: start state, 2: end state, 3: total
[lossTransitionProbability2D,lossTransitionProbabilityMats] = GET_TRANS_PROBS(task_1, subjInd, start);
[lossTransitionProbabilityNoPersist2D,lossTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(task_1, subjInd, start);
                                                                              % calculate transition probabilities while excluding state persistence, i.e. independent of autocorrelation
[rewardTransitionProbability2D,rewardTransitionProbabilityMats] = GET_TRANS_PROBS(task_2, subjInd, start);
[rewardTransitionProbabilityNoPersist2D,rewardTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(task_2, subjInd, start);

%
partitions_con = reshape(partitions_tmp, [1,numel(partitions_tmp)]);

subjInd_all = [];
for si = 1:Nsubj
    subjInd_all = [subjInd_all; ones(size(partitions_tmp,2), 1)*si];
end

[wholeTransitionProbability2D,wholeTransitionProbabilityMats] = GET_TRANS_PROBS(partitions_con, subjInd_all, start);
[wholeTransitionProbabilityNoPersist2D,wholeTransitionProbabilityMatsNoPersist] = GET_TRANS_PROBS_NO_PERSIST(partitions_con, subjInd_all, start);

label_idx = cell(K,1);
for kk = linspace(1,K,K)
    label_idx{kk} = num2str(kk);
end

figure; imagesc(squeeze(mean(wholeTransitionProbabilityMatsNoPersist,1)) .* ~eye(K));
xticks(1:K); yticks(1:K); colormap('plasma'); 
xticklabels(label_idx); yticklabels(label_idx); axis square;
COLOR_TICK_LABELS(true,true,K);
ylabel('Current State'); xlabel('Next New State');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
cmap = fcn_cmaphot; colormap(cmap);colorbar;
set(gcf,'color','w');

figure; imagesc(squeeze(mean(wholeTransitionProbabilityMats,1)));
xticks(1:K); yticks(1:K); colormap('plasma'); 
xticklabels(label_idx); yticklabels(label_idx); axis square;
COLOR_TICK_LABELS(true,true,K);
ylabel('Current State'); xlabel('Next New State');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
cmap = fcn_cmaphot; colormap(cmap);colorbar;
set(gcf,'color','w');
%% permutation testing to compare transition probability matrices
disp('start permutation testing')
nperms = 100000;
[~,pvalrew_lo] = PERM_TEST(rewardTransitionProbabilityMatsNoPersist,lossTransitionProbabilityMatsNoPersist,nperms);
[~,pvallo_rew] = PERM_TEST(lossTransitionProbabilityMatsNoPersist,rewardTransitionProbabilityMatsNoPersist,nperms);

% plot transition probabilities
grpAvgLoss = squeeze(mean(lossTransitionProbabilityMatsNoPersist,1)) .* ~eye(K);
% nans occur for transitions from states that are not present at all for a subject or within the tested blocks. 
% transitions to that state are 0
grpAvgReward = squeeze(nanmean(rewardTransitionProbabilityMatsNoPersist,1)) .* ~eye(K); 

maxVal = max(max([grpAvgLoss,grpAvgReward])); % sync color scales
%%
label_idx = cell(K,1);
for kk = linspace(1,K,K)
    label_idx{kk} = num2str(kk);
end

f1 = figure;

subplot(2,2,1);
imagesc(grpAvgLoss);
xticks(1:K); yticks(1:K); colormap('plasma'); 
xticklabels(label_idx); yticklabels(label_idx); axis square;
COLOR_TICK_LABELS(true,true,K);
ylabel('Current State'); xlabel('Next New State');
title('Loss');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colormap(cmap);colorbar;

%
subplot(2,2,2);
imagesc(grpAvgReward);
xticks(1:K); yticks(1:K); colormap('plasma'); 
xticklabels(label_idx); yticklabels(label_idx); axis square;
COLOR_TICK_LABELS(true,true,K);
ylabel('Current State'); xlabel('Next New State');
title('Reward');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colormap(cmap);colorbar;

%
subplot(2,2,3);
RewardMinusLossTP = (grpAvgReward-grpAvgLoss);
imagesc(RewardMinusLossTP); colormap('plasma');
xticks(1:K); xticklabels(label_idx); 
yticks(1:K); yticklabels(label_idx); axis square
ylabel('Current State'); xlabel('Next New State');
sig_thresh = 0.05 / K^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvalrew_lo < sig_thresh);
for si = 1:length(x)
    if y(si) ~= x(si)
        text(x(si)-.12,y(si)+.12,'*','Color','b');
    end
end
% caxis_bound = max(max(abs(RewardMinusLossTP)));
% h = colorbar; 
% ylabel(h,'reward-loss'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,K);
title('Reward > Loss');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
colormap(cmap);colorbar;

%
subplot(2,2,4);
LossMinusRewardTP = (grpAvgLoss-grpAvgReward);
imagesc(LossMinusRewardTP); colormap('plasma');
xticks(1:K); xticklabels(label_idx); 
yticks(1:K); yticklabels(label_idx); axis square
ylabel('Current State'); xlabel('Next New State');
sig_thresh = 0.05 / K^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvallo_rew < sig_thresh);
for si = 1:length(x)
    if y(si) ~= x(si)
        text(x(si)-.12,y(si)+.12,'*','Color','b');
    end
end
%caxis_bound = max(max(abs(LossMinusRewardTP)));
%h = colorbar; 
%ylabel(h,'loss-reward'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,K);
title('Loss > Reward');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
colormap(cmap);colorbar;
set(gcf,'color','w');

%%
DT_lo = LT1;
DT_re = LT2;
FO_lo = P1;
FO_re = P2;
TP_lo = lossTransitionProbabilityNoPersist2D;
TP_re = rewardTransitionProbabilityNoPersist2D;

cd([root_path, 'results/'])
check = exist('graph_met', 'dir')
if check == 0
    mkdir('graph_met')
end
save([root_path, 'results/graph_met/DT_FO_TP.mat'], 'DT_lo', 'DT_re', 'FO_lo', 'FO_re', 'TP_lo', 'TP_re')