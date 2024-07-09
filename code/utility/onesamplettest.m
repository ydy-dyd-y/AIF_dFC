% one sample ttest
clear
clc
close all

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

mi = 4; str = sprintf('%02d', mi); % model index

% load 32k surfaces
addpath(genpath([root_path,'tools/']));
load([root_path,'tools/fcn/surfinfo.mat'])
n_net = 7;
aa = niftiread([root_path,'tools/fcn/Schaefer2018_100Parcels_',num2str(n_net),'Networks_order.dlabel.nii']);
aa = squeeze(aa);
ar = aa(32493:end); al = aa(1:32492);
cr = zeros(size(ar)); cl = zeros(size(al));
cmap = fcn_cmaphot;

region_level = 2; % 1:node, 2:edge, 3:state
behav_data_type = 1; % 1: original hrf % dISC of hrf
thres = 0.01;
pq = 'p';

plot_method = 2; % one color, 2: two color
region_level_label = {'node', 'edge', 'state'};
method_label = {'act', 'dISC'};
data_type = [method_label{behav_data_type}, '_', region_level_label{region_level}, '_DS'];
results_path = [root_path, 'results/GLM']; cd(results_path)

load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);
originLabels = {'R-HIP', 'R-AMY', 'R-pTHA', 'R-aTHA', 'R-NAc', 'R-PUT', 'R-CAU', ...
    'L-HIP', 'L-AMY', 'L-pTHA', 'L-aTHA', 'L-NAc', 'L-PUT', 'L-CAU', ...
    'L-GPe', 'L-GPi', 'L-SNc', 'L-red', 'L-SNr', 'L-PBPN', 'L-VTA', 'L-VP', 'L-HABN', 'L-hypoTHAL',  'L-MM', 'L-STN', ...
    'R-GPe', 'R-GPi', 'R-SNc', 'R-red', 'R-SNr', 'R-PBPN', 'R-VTA', 'R-VP', 'R-HABN', 'R-hypoTHAL', 'R-MM', 'R-STN', ...
    'DR', 'LC', 'LDTg', 'MnR', 'mRt', 'PAG', 'PBC', 'PnO', 'PTg'};
wholeLabels = [cortex_label; originLabels'];

load([root_path, 'data/connectivity_sign.mat'], 'H')

Nnodes = numel(wholeLabels);
Nedges = Nnodes*(Nnodes - 1)/2;

ir = [];
ic = [];
for i = 1:Nnodes
    ir = [ir;repelem(i, Nnodes - i)'];
    ic = [ic;linspace(i+1, Nnodes, Nnodes - i)'];
end

H_vec = [];
for ii = 1:length(ir)
    H_vec = [H_vec; H(ir(ii), ic(ii))];
end

for ii = 1:Nedges
    ROE_labels{1,ii} = wholeLabels{ir(ii)};
    ROE_labels{2,ii} = wholeLabels{ic(ii)};
end

if n_net == 7
    net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SC'};
elseif n_net == 17
    net = {'VISCent', 'VISPer', 'SOMa','SOMb', 'DATa', 'DATb', 'VATa', 'VATb', 'LIMa', 'LIMb', 'ContA', 'ContB', 'ContC', 'DMNa', 'DMNb', 'DMNc', 'TP', 'SC'};
end

labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
lab = labs(1:100);
lab_sub = [(2:8)'; (2:8)'; (10:2:32)'; (9:2:31)'; ones(9,1)*33] + 6;
lab = [lab; lab_sub];
[gx,gy,idx] = grid_communities(lab); % BCT function
subLabels = {'R-HIP','L-HIP','R-AMY','L-AMY','R-pTHA','L-pTHA','R-aTHA','L-aTHA', 'R-NAc','L-NAc','R-PUT','L-PUT','R-CAU','L-CAU',...
    'R-GPe','L-GPe','R-GPi','L-GPi','R-SNc','L-SNc','R-red','L-red','R-SNr','L-SNr','R-PBPN','L-PBPN',...
    'R-VTA','L-VTA', 'R-VP','L-VP','R-HABN','L-HABN','R-hypoTHAL','L-hypoTHAL','R-MM','L-MM', 'R-STN','L-STN', ...
    'AAN'};

tmp_name = [data_type, num2str(str), '.mat'];

load(tmp_name, 'params_mass')   % (147    54     EV ) : (Nnodes, Nsubj, EV)
Nregions = size(params_mass, 1);
Nsubj = size(params_mass, 2);
num_var = size(params_mass, 3);

mean_observed = squeeze(mean(params_mass, 2));% (147 x EV x 500)
%%
p_rt = []; sig_p_rt = [];
p_lt = []; sig_p_lt = [];
for roi = 1:Nregions
    for vi = 1:num_var
        [sig_p_rt(roi,vi), p_rt(roi,vi)] = ttest(params_mass(roi, :, vi), 0, 'Alpha', thres, 'Tail', 'right');
        [sig_p_lt(roi,vi), p_lt(roi,vi)] = ttest(params_mass(roi, :, vi), 0, 'Alpha', thres, 'Tail', 'left');
    end
end
q_rt = [];
q_lt = [];
for vi = 1:num_var
    q_tmp = mafdr(p_rt(:,vi),'BHFDR', true);
    q_rt = [q_rt, q_tmp];
    q_tmp = mafdr(p_lt(:,vi),'BHFDR', true);
    q_lt = [q_lt, q_tmp];
end
sig_q_rt = q_rt < thres;
sig_q_lt = q_lt < thres;

Label_sig_p_exc_pos = {};
Label_sig_p_inh_pos = {};
Label_sig_p_exc_neg = {};
Label_sig_p_inh_neg = {};
Label_sig_q_exc_pos = {};
Label_sig_q_inh_pos = {};
Label_sig_q_exc_neg = {};
Label_sig_q_inh_neg = {};
for vi = 1:num_var
    i_ep = 1; i_ip = 1; i_en = 1; i_in = 1;
    for roi = 1:Nregions
        if sig_p_rt(roi, vi) == 1
            if H_vec(roi) == 1
                Label_sig_p_exc_pos{i_ep,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_ep = i_ep+1;
            elseif H_vec(roi) == -1
                Label_sig_p_inh_pos{i_ip,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_ip = i_ip+1;
            end
        end
        if sig_p_lt(roi, vi) == 1
            if H_vec(roi) == 1
                Label_sig_p_exc_neg{i_en,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_en = i_en+1;
            elseif H_vec(roi) == -1
                Label_sig_p_inh_neg{i_in,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_in = i_in+1;
            end
        end
    end
    i_ep = 1; i_ip = 1; i_en = 1; i_in = 1;
    for roi = 1:Nregions
        if sig_q_rt(roi, vi) == 1
            if H_vec(roi) == 1
                Label_sig_q_exc_pos{i_ep,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_ep = i_ep+1;
            elseif H_vec(roi) == -1
                Label_sig_q_inh_pos{i_ip,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_ip = i_ip+1;
            end
        end
        if sig_q_lt(roi, vi) == 1
            if H_vec(roi) == 1
                Label_sig_q_exc_neg{i_en,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_en = i_en+1;
            elseif H_vec(roi) == -1
                Label_sig_q_inh_neg{i_in,vi} = [ROE_labels{1, roi}, ' : ', ROE_labels{2, roi}];
                i_in = i_in+1;
            end
        end
    end
end

%%
load([root_path, 'data/cm.mat'], 'cm')
if mi == 1
    var_idx = 6;
elseif ismember(mi, [4,5,6,7])
    var_idx = 5;
end
beta_map_exc_pos = zeros(Nnodes);
beta_map_exc_neg = zeros(Nnodes);
beta_map_inh_pos = zeros(Nnodes);
beta_map_inh_neg = zeros(Nnodes);

beta_exc_pos = zeros(Nnodes, Nnodes, Nsubj);
beta_exc_neg = zeros(Nnodes, Nnodes, Nsubj);
beta_inh_pos = zeros(Nnodes, Nnodes, Nsubj);
beta_inh_neg = zeros(Nnodes, Nnodes, Nsubj);

beta_exc = zeros(Nnodes, Nnodes, Nsubj);
beta_inh = zeros(Nnodes, Nnodes, Nsubj);

if pq == 'p'
    for ii = 1:Nedges
        if H_vec(ii) == 1
            beta_map_exc_pos(ir(ii), ic(ii)) = sig_p_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_pos(ic(ii), ir(ii)) = sig_p_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_neg(ir(ii), ic(ii)) = sig_p_lt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_neg(ic(ii), ir(ii)) = sig_p_lt(ii, var_idx) * mean_observed(ii, var_idx);
            for si = 1:Nsubj
                if params_mass(ii, si, var_idx) > 0
                    beta_exc_pos(ir(ii), ic(ii), si) = sig_p_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc_pos(ic(ii), ir(ii), si) = sig_p_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc(ir(ii), ic(ii), si) = params_mass(ii, si, var_idx);
                    beta_exc(ic(ii), ir(ii), si) = params_mass(ii, si, var_idx);
                end
                if params_mass(ii, si, var_idx) < 0
                    beta_exc_neg(ir(ii), ic(ii), si) = -sig_p_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc_neg(ic(ii), ir(ii), si) = -sig_p_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                end
                
            end
        elseif H_vec(ii) == -1
            beta_map_inh_pos(ir(ii), ic(ii)) = sig_p_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_pos(ic(ii), ir(ii)) = sig_p_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_neg(ir(ii), ic(ii)) = sig_p_lt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_neg(ic(ii), ir(ii)) = sig_p_lt(ii, var_idx) * mean_observed(ii, var_idx);
            for si = 1:Nsubj
                if params_mass(ii, si, var_idx) > 0
                    beta_inh_pos(ir(ii), ic(ii), si) = sig_p_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh_pos(ic(ii), ir(ii), si) = sig_p_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                end
                if params_mass(ii, si, var_idx) < 0
                    beta_inh_neg(ir(ii), ic(ii), si) = -sig_p_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh_neg(ic(ii), ir(ii), si) = -sig_p_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh(ir(ii), ic(ii), si) = -params_mass(ii, si, var_idx);
                    beta_inh(ic(ii), ir(ii), si) = -params_mass(ii, si, var_idx);
                end
                
            end
        end
    end
elseif pq == 'q'
    for ii = 1:Nedges
        if H_vec(ii) == 1
            beta_map_exc_pos(ir(ii), ic(ii)) = sig_q_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_pos(ic(ii), ir(ii)) = sig_q_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_neg(ir(ii), ic(ii)) = sig_q_lt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_exc_neg(ic(ii), ir(ii)) = sig_q_lt(ii, var_idx) * mean_observed(ii, var_idx);
            for si = 1:Nsubj
                if params_mass(ii, si, var_idx) > 0
                    beta_exc_pos(ir(ii), ic(ii), si) = sig_q_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc_pos(ic(ii), ir(ii), si) = sig_q_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc(ir(ii), ic(ii), si) = params_mass(ii, si, var_idx);
                    beta_exc(ic(ii), ir(ii), si) = params_mass(ii, si, var_idx);
                end
                if params_mass(ii, si, var_idx) < 0
                    beta_exc_neg(ir(ii), ic(ii), si) = -sig_q_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_exc_neg(ic(ii), ir(ii), si) = -sig_q_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                end
                
            end
        elseif H_vec(ii) == -1
            beta_map_inh_pos(ir(ii), ic(ii)) = sig_q_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_pos(ic(ii), ir(ii)) = sig_q_rt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_neg(ir(ii), ic(ii)) = sig_q_lt(ii, var_idx) * mean_observed(ii, var_idx);
            beta_map_inh_neg(ic(ii), ir(ii)) = sig_q_lt(ii, var_idx) * mean_observed(ii, var_idx);
            for si = 1:Nsubj
                if params_mass(ii, si, var_idx) > 0
                    beta_inh_pos(ir(ii), ic(ii), si) = sig_q_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh_pos(ic(ii), ir(ii), si) = sig_q_rt(ii, var_idx) * params_mass(ii, si, var_idx);
                end
                if params_mass(ii, si, var_idx) < 0
                    beta_inh_neg(ir(ii), ic(ii), si) = -sig_q_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh_neg(ic(ii), ir(ii), si) = -sig_q_lt(ii, var_idx) * params_mass(ii, si, var_idx);
                    beta_inh(ir(ii), ic(ii), si) = -params_mass(ii, si, var_idx);
                    beta_inh(ic(ii), ir(ii), si) = -params_mass(ii, si, var_idx);
                end
            end
        end
    end
end

figure;
ax(1) = subplot(1,4,1);
imagesc(beta_map_exc_pos(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
colorbar('southoutside');
if n_net == 7
    for i = 1:numel(net)-1
        labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
        labelhandle.HorizontalAlignment = 'center';
    end
    
    for i = 1:numel(subLabels)
        if subLabels{i}(1:2) == 'R-'
            sublabelhandle = text(Nnodes+2, 100.5+i, extractAfter(subLabels{i}, 'R-'));
            set(sublabelhandle, 'FontSize', 5.2)
        elseif subLabels{i}(1:2) == 'L-'
            continue
        else
            sublabelhandle = text(Nnodes+2, 100+i, subLabels{i});
            set(sublabelhandle, 'FontSize', 4.6)
        end
    end
end
ax(2) = subplot(1,4,2);
imagesc(beta_map_exc_neg(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
colorbar('southoutside');
if n_net == 7
    for i = 1:numel(net)-1
        labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
        labelhandle.HorizontalAlignment = 'center';
    end
    
    for i = 1:numel(subLabels)
        if subLabels{i}(1:2) == 'R-'
            sublabelhandle = text(Nnodes+2, 100.5+i, extractAfter(subLabels{i}, 'R-'));
            set(sublabelhandle, 'FontSize', 5.2)
        elseif subLabels{i}(1:2) == 'L-'
            continue
        else
            sublabelhandle = text(Nnodes+2, 100+i, subLabels{i});
            set(sublabelhandle, 'FontSize', 4.6)
        end
    end
end
ax(3) = subplot(1,4,3);
imagesc(beta_map_inh_pos(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
colorbar('southoutside');
if n_net == 7
    for i = 1:numel(net)-1
        labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
        labelhandle.HorizontalAlignment = 'center';
    end
    for i = 1:numel(subLabels)
        if subLabels{i}(1:2) == 'R-'
            sublabelhandle = text(Nnodes+2, 100.5+i, extractAfter(subLabels{i}, 'R-'));
            set(sublabelhandle, 'FontSize', 5.2)
        elseif subLabels{i}(1:2) == 'L-'
            continue
        else
            sublabelhandle = text(Nnodes+2, 100+i, subLabels{i});
            set(sublabelhandle, 'FontSize', 4.6)
        end
    end
end
ax(4) = subplot(1,4,4);
imagesc(beta_map_inh_neg(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
colorbar('southoutside');
if n_net == 7
    for i = 1:numel(net)-1
        labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
        labelhandle.HorizontalAlignment = 'center';
    end
    
    for i = 1:numel(subLabels)
        if subLabels{i}(1:2) == 'R-'
            sublabelhandle = text(Nnodes+2, 100.5+i, extractAfter(subLabels{i}, 'R-'));
            set(sublabelhandle, 'FontSize', 5.2)
        elseif subLabels{i}(1:2) == 'L-'
            continue
        else
            sublabelhandle = text(Nnodes+2, 100+i, subLabels{i});
            set(sublabelhandle, 'FontSize', 4.6)
        end
    end
end

colormap(ax(1), cm(102:end,:))
colormap(ax(2), cm(1:101,:))
colormap(ax(3), cm(102:end,:))
colormap(ax(4), cm(1:101,:))
set(gcf, 'color', 'w')

if region_level == 2
    inter_whole = cell({});
    con_eff = cell({'exc_pos', 'inh_pos', 'exc_neg', 'inh_neg'});
    for ci = 1:length(con_eff)
        eval(['beta_tmp = beta_', con_eff{ci}, ';'])
        inter = [];
        for si = 1:Nsubj
            [inter_tmp, intra_tmp] = integration_recruitment(squeeze(beta_tmp(:, :, si)), labs);
            inter_tmp = inter_tmp - diag(diag(inter_tmp));
            inter_tmp = inter_tmp + diag(diag(intra_tmp));
            inter = cat(3, inter, inter_tmp);
        end
        inter_whole{ci} = inter;
    end
    inter_whole_mean = cellfun(@(x) nanmean(x, 3), inter_whole, 'UniformOutput', false);
    
    check = exist([root_path, 'results/graph_met/'], 'dir');
    if check == 0
        mkdir([root_path, 'results/graph_met/'])
    end
    save([root_path, 'results/graph_met/bct_metric_per4contype_DS', sprintf('%02d', mi), '.mat'], 'inter_whole', 'inter_whole_mean')
    
    inter_whole = cell({});
    beta_tmp = beta_inh + beta_exc;
    inter = [];
    for si = 1:Nsubj
        [inter_tmp, intra_tmp] = integration_recruitment(squeeze(beta_tmp(:, :, si)), labs);
        inter_tmp = inter_tmp - diag(diag(inter_tmp));
        inter_tmp = inter_tmp + diag(diag(intra_tmp));
        inter = cat(3, inter, inter_tmp);
    end
    inter_whole{1} = inter;
    inter_whole_mean = cellfun(@(x) nanmean(x, 3), inter_whole, 'UniformOutput', false);
    
    check = exist([root_path, 'results/graph_met/'], 'dir');
    if check == 0
        mkdir([root_path, 'results/graph_met/'])
    end
    save([root_path, 'results/graph_met/bct_metric_DS', sprintf('%02d', mi), '.mat'], 'inter_whole', 'inter_whole_mean')

end