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
%%
kgL = []; kgR = [];
peL = []; peR = [];

for i = 1:Nsubj
    
    est_kf_lore{i,1} = tapas_fitModel(data{i}.response{1}, data{i}.input{1}, 'tapas_kf_config','tapas_softmax_binary_config');
    est_kf_lore{i,2} = tapas_fitModel(data{i}.response{2}, data{i}.input{2}, 'tapas_kf_config','tapas_softmax_binary_config');

    kgL = [kgL, est_kf_lore{i,1}.traj.g];
    kgR = [kgR, est_kf_lore{i,2}.traj.g]; 
    
    peL = [peL, est_kf_lore{i,1}.traj.da];
    peR = [peR, est_kf_lore{i,2}.traj.da]; 
    
end

save([root_path, 'results/behav/kalman.mat'], 'est_kf_lore')
rm_sub = [22, 36, 44, 57];

%%
plot(mean(kgL,2)); 
hold on
plot(mean(kgR,2));

legend('Loss', 'Reward');

figure; heatmap(peL')

%%
idx_sub = 5;
tapas_kf_plotTraj(est_kf_lore{idx_sub, 1});
figure; plot(kgL(:,idx_sub))
