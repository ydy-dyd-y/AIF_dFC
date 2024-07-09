%% Plotting Supplementary Fig. 1
clc
clear

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

addpath([root_path,'tools/']);
behv_path = [root_path, 'data/behav/raw_behav_data/'];
data_list = {dir([behv_path, 'YAD*.csv']).name};
Nsubj = numel(data_list);
Nsubidx = str2double(data_list{end}(7:9));

for i = 1:Nsubidx
    sname{i} = ['YAD_1', sprintf('%04d', i)];
end

data = preproc_behave(behv_path,sname,1:Nsubidx);
tmp = cellfun(@(i) isempty(i), data, 'UniformOutput', false); tmp = cell2mat(tmp);
data(find(tmp)) = [];
save([root_path, 'data/behav/behavior_data.mat'], 'data')
%%
% loss_resp : loss task에서 (trial X subjects)로 응답한 답 matrix
% reward_resp : reward task에서 (trial X subjects)로 응답한 답 matrix

loss_resp = [];
reward_resp = [];
for i = 1:Nsubj
    if isempty(data{i})
        continue
    end
    loss_resp = [loss_resp, data{i,1}.response{1,1}];
    reward_resp = [reward_resp, data{i,1}.response{1,2}];
end     
%% loss_corr, reward_corr 각각, 행렬의 요소가 맞았는지 틀렸는지로 바뀜(틀린건 0, 맞은건 1)
loss_corr = [];
reward_corr = [];
for i = 1:numel(data_list)
    for k = 1:size(loss_resp, 1)
        if loss_resp(k,i) == data{4, 1}.input{1,1}(k)
            loss_corr(k,i) = 1;
        else loss_corr(k,i) = 0 ;
        end
    end
end

for i = 1:numel(data_list)
    for k = 1:size(reward_resp, 1)
        if reward_resp(k,i) == data{4, 1}.input{1,2}(k)
            reward_corr(k,i) = 1;
        else reward_corr(k,i) = 0 ;
        end
    end
end 

%% wrong_1, wrong_2는 각각 틀린 문제의 번호를 모든 사람에서 모은 array
wrong_1 = [];
wrong_2 = [];

for i = 1:numel(data_list)
   a = find(loss_corr(:,i) == 0);
   b = find(reward_corr(:,i) == 0);
   wrong_1 = cat(1, wrong_1, a);
   wrong_2 = cat(1, wrong_2, b);
   
end

%% find wrong trials index for each subject
wrong_lo = zeros(numel(data_list), size(loss_resp, 1));
wrong_re = zeros(numel(data_list), size(reward_resp, 1));

for i = 1:numel(data_list)
   a = find(loss_corr(:,i) == 0);
   b = find(reward_corr(:,i) == 0);
   wrong_lo(i,1:length(a)) = a;
   wrong_re(i,1:length(b)) = b;
end

%%
f = figure;
subplot(1,2,1)
histogram(wrong_1);
title('Loss task')
xlabel('Trial')
ylabel('Number of subjects who lose')

subplot(1,2,2)
histogram(wrong_2);
title('Reward task');
xlabel('Trial')

set(gcf, 'color', 'w')
set(gcf, 'position', [840 439 930 296])
savefile = fullfile([root_path, 'results/Figures/'], 'supp_fig1.pdf');
exportgraphics(f, savefile)