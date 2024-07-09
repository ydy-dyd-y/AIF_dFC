% fig. 7a~c with (mi = 13, behav_data_type = 1, region_level = 2)
clear
clc
close all

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

mi = 8; str = sprintf('%02d', mi); % model index

region_level = 2; behav_data_type = 1;
region_level_label = {'node', 'edge', 'state', 'system'};
method_label = {'act', 'dISC'};
data_type = [method_label{behav_data_type}, '_', region_level_label{region_level}, '_DS'];

n_net = 7;
load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);
originLabels = {'R-HIP', 'R-AMY', 'R-pTHA', 'R-aTHA', 'R-NAc', 'R-PUT', 'R-CAU', ...
    'L-HIP', 'L-AMY', 'L-pTHA', 'L-aTHA', 'L-NAc', 'L-PUT', 'L-CAU', ...
    'L-GPe', 'L-GPi', 'L-SNc', 'L-red', 'L-SNr', 'L-PBPN', 'L-VTA', 'L-VP', 'L-HABN', 'L-hypoTHAL',  'L-MM', 'L-STN', ...
    'R-GPe', 'R-GPi', 'R-SNc', 'R-red', 'R-SNr', 'R-PBPN', 'R-VTA', 'R-VP', 'R-HABN', 'R-hypoTHAL', 'R-MM', 'R-STN', ...
    'DR', 'LC', 'LDTg', 'MnR', 'mRt', 'PAG', 'PBC', 'PnO', 'PTg'};
wholeLabels = [cortex_label; originLabels'];

% load 32k surfaces
addpath(genpath([root_path,'tools/']));
load([root_path,'tools/fcn/surfinfo.mat'])
aa = niftiread([root_path,'tools/fcn/Schaefer2018_100Parcels_',num2str(n_net),'Networks_order.dlabel.nii']);
load([root_path, 'data/cm.mat'], 'cm')

Nnodes = numel(wholeLabels);
Nedges = Nnodes * (Nnodes - 1)/2;
ir = [];
ic = [];
for i = 1:Nnodes
    ir = [ir;repelem(i, Nnodes - i)'];
    ic = [ic;linspace(i+1, Nnodes, Nnodes - i)'];
end
if n_net == 7
    net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SC'};
    if region_level == 5
        net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN_1', 'DMN_2', 'DMN_3'};
    end
elseif n_net == 17
    net = {'VISCent', 'VISPer', 'SOMa','SOMb', 'DATa', 'DATb', 'VATa', 'VATb', 'LIMa', 'LIMb', 'ContA', 'ContB', 'ContC', 'DMNa', 'DMNb', 'DMNc', 'TP', 'SC'};
end
labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
lab = labs(1:100);
%lab(101:end) = [(2:8)'; (2:8)'; (9:20)'; (9:20)'; (21:29)'] + (n_net - 1);   % length of gx increase from 48 -> 210
[gx,gy,idx] = grid_communities(lab); % BCT function
subLabels = {'R-HIP','L-HIP','R-AMY','L-AMY','R-pTHA','L-pTHA','R-aTHA','L-aTHA', 'R-NAc','L-NAc','R-PUT','L-PUT','R-CAU','L-CAU',...
    'R-GPe','L-GPe','R-GPi','L-GPi','R-SNc','L-SNc','R-red','L-red','R-SNr','L-SNr','R-PBPN','L-PBPN',...
    'R-VTA','L-VTA', 'R-VP','L-VP','R-HABN','L-HABN','R-hypoTHAL','L-hypoTHAL','R-MM','L-MM', 'R-STN','L-STN', ...
    'AAN'};
gx2 = [46.5, 46.5, 100.5, 100.5, 46.5];
gy2 = [46.5, 100.5, 100.5, 46.5, 46.5];

cd([root_path, 'results/GLM'])
tmp_name = [data_type, num2str(str), '.mat'];
if ismember(mi, [1])
    var_position = 4:5;
elseif ismember(mi, [4])
    var_position = 3:6;
elseif ismember(mi, [2, 3, 5])
    var_position = 3:4;
elseif ismember(mi, [7, 8, 9])
    var_position = 3:9
elseif mi == 6
    var_position = 3:5;
end

load(tmp_name, 'params_mass', 'p_values_mass')
Nsubj = size(params_mass, 2);
mean_observed = squeeze(mean(params_mass,2));
X = [];
for var_positioni = var_position
    X = cat(3, X, params_mass(:,:,var_positioni));
end

Nregions = size(X, 1);
%% defining parameters
bc = 0; pq = 'p';
if bc == 1
    thres = 0.05/Nregions;
else
    thres = 0.01;
end

if length(var_position) == 3
    x1_idx_cell = {[1,2], [1,3], [2,3]};
    ft_name = {'certain', 'uncertain1', 'uncertain2'};
elseif length(var_position) == 4
    x1_idx_cell = {[1,2], [2,3], [3,4], [1,3], [1,4], [2,4]};
    ft_name = {'pw', 'pl', 'nw', 'nl'};
elseif length(var_position) == 7
    x1_idx_cell = {[1,2], [1,3], [1,4], [1,5], [1,6], [1,7]};
    ft_name = {'certain', 'uc1', 'uc2_1', 'uc2_2', 'uc2_3', 'uc2_4', 'uc2_5'};
elseif length(var_position) == 2
    x1_idx_cell = {[1,2]};
    if mi == 1
        ft_name = {'AS', 'PREC'};
    else
        ft_name = {'Certain', 'Uncertain'};
    end
end

if region_level == 1
    Label_sig_rt = {}; Label_sig_lt = {};
end

iii = 1;
for compar = 1:length(x1_idx_cell)
    if region_level == 2
        eval(['f_whole_set', num2str(compar), '= figure;'])
        ax_i = 1;
    end
    Sig_rt_mat = []; Sig_rt_cor_mat =[];
    Sig_lt_mat = []; Sig_lt_cor_mat =[];
    for x1_idx = x1_idx_cell{compar}
        x2_idx = setdiff(x1_idx_cell{compar}, x1_idx);
        h_rt = []; p_rt = []; h_rt_cor = []; p_rt_cor = [];
        h_lt = []; p_lt = []; h_lt_cor = []; p_lt_cor = [];
        %H_rt = []; P_rt = []; H_lt = []; P_lt = [];
        for ni = 1:Nregions
            if region_level == 2
                if mean_observed(ni,var_position(x1_idx)) > 0
                    [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', thres, 'Tail', 'right');  % right : x1 > x2
                    h_rt = [h_rt, h]; p_rt = [p_rt, p];
                else
                    h = 0; p = 1; h_rt = [h_rt, h]; p_rt = [p_rt, p];
                end
                if mean_observed(ni,var_position(x1_idx)) < 0
                    [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', thres, 'Tail', 'left');  % left : x1 < x2
                    h_lt = [h_lt, h]; p_lt = [p_lt, p];
                else
                    h = 0; p = 1; h_lt = [h_lt, h]; p_lt = [p_lt, p];
                end
            else
                
                [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', thres, 'Tail', 'right');  % right : x1 > x2
                h_rt = [h_rt, h]; p_rt = [p_rt, p];
                [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', 0.05/(n_net*(n_net+1)/2), 'Tail', 'right');
                h_rt_cor = [h_rt_cor, h]; p_rt_cor = [p_rt_cor, p];
                
                [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', thres, 'Tail', 'left');   % left : x1 < x2
                h_lt = [h_lt, h]; p_lt = [p_lt, p];
                [h, p] = ttest(X(ni, :, x1_idx), X(ni, :, x2_idx), 'Alpha', 0.05/(n_net*(n_net+1)/2), 'Tail', 'left');
                h_lt_cor = [h_lt_cor, h]; p_lt_cor = [p_lt_cor, p];
                
            end
        end
        P_rt{iii} = p_rt; H_rt{iii} = h_rt; H_rt_cor{iii} = h_rt_cor;
        P_lt{iii} = p_lt; H_lt{iii} = h_lt; H_lt_cor{iii} = h_lt_cor;
        
        if region_level == 1
            p_sig_node_rt{iii} = find(h_rt); p_sig_node_lt{iii} = find(h_lt);
            Label_sig_rt = {}; Label_sig_lt = {};
            if pq == 'p'
                mean_beta_sig_rt = mean_observed(p_sig_node_rt, var_position(x1_idx));
                mean_beta_sig_lt = mean_observed(p_sig_node_lt, var_position(x1_idx));
                for roi = 1:length(p_sig_node_rt)
                    Label_sig_rt{1,roi} = wholeLabels{ir(p_sig_node_rt(roi))};
                    Label_sig_rt{2,roi} = wholeLabels{ic(p_sig_node_rt(roi))};
                end
                for roi = 1:length(p_sig_node_lt)
                    Label_sig_lt{1,roi} = wholeLabels{ir(p_sig_node_lt(roi))};
                    Label_sig_lt{2,roi} = wholeLabels{ic(p_sig_node_lt(roi))};
                end
            elseif pq == 'q'
                mean_beta_sig_rt = mean_observed(q_sig_node_rt, var_position(x1_idx));
                mean_beta_sig_lt = mean_observed(q_sig_node_lt, var_position(x1_idx));
                for roi = 1:length(q_sig_node_rt)
                    Label_sig_rt{1,roi} = wholeLabels{ir(q_sig_node_rt(roi))};
                    Label_sig_rt{2,roi} = wholeLabels{ic(q_sig_node_rt(roi))};
                end
                for roi = 1:length(q_sig_node_lt)
                    Label_sig_lt{1,roi} = wholeLabels{ir(q_sig_node_lt(roi))};
                    Label_sig_lt{2,roi} = wholeLabels{ic(q_sig_node_lt(roi))};
                end
            end
        elseif region_level == 2
            beta_map_rt = zeros(Nnodes);
            beta_map_lt = zeros(Nnodes);
            if pq == 'p'
                for ii = 1:Nedges
                    beta_map_rt(ir(ii), ic(ii)) = h_rt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_rt(ic(ii), ir(ii)) = h_rt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_lt(ir(ii), ic(ii)) = h_lt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_lt(ic(ii), ir(ii)) = h_lt(ii) * mean_observed(ii, var_position(x1_idx));
                end
            elseif pq == 'q'
                for ii = 1:Nedges
                    beta_map_rt(ir(ii), ic(ii)) = h_rt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_rt(ic(ii), ir(ii)) = h_rt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_lt(ir(ii), ic(ii)) = h_lt(ii) * mean_observed(ii, var_position(x1_idx));
                    beta_map_lt(ic(ii), ir(ii)) = h_lt(ii) * mean_observed(ii, var_position(x1_idx));
                end
            end
            eval(['figure(f_whole_set', num2str(compar), ')'])
            ax(ax_i) = subplot(2, 2, ax_i);
            ax_i = ax_i + 1;
            imagesc(beta_map_rt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
            colorbar('southoutside');
            if n_net == 7
                for i = 1:numel(net)-1
                    labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
                    labelhandle.HorizontalAlignment = 'center';
                end
            end
            title([ft_name{x1_idx}, ' > ', ft_name{x2_idx}])
            ax(ax_i) = subplot(2, 2, ax_i);
            ax_i = ax_i + 1;
            imagesc(beta_map_lt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
            colorbar('southoutside');
            if n_net == 7
                for i = 1:numel(net)-1
                    labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
                    labelhandle.HorizontalAlignment = 'center';
                end
            end
            title([ft_name{x1_idx}, ' < ', ft_name{x2_idx}])
            colormap(ax(ax_i-2), cm(102:end,:))
            colormap(ax(ax_i-1), cm(1:101,:))
            set(gcf, 'color', 'w')
            
            fig_name = ['set_',num2str(compar),'f_save', num2str(ax_i - 2)];
            eval([fig_name, ' = figure;'])
            imagesc(beta_map_rt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
            colorbar('southoutside'); colormap(cm(102:end,:))
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
            set(gcf, 'color', 'w')
            
            fig_name = ['set_',num2str(compar),'f_save', num2str(ax_i - 1)];
            eval([fig_name, ' = figure;'])
            imagesc(beta_map_lt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5);
            colorbar('southoutside'); colormap(cm(1:101,:))
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
            set(gcf, 'color', 'w')
            
            beta_tmp = zeros(Nnodes);
            for ii = 1:Nedges
                beta_tmp(ir(ii), ic(ii)) = h_rt(ii);
                beta_tmp(ic(ii), ir(ii)) = h_rt(ii);
            end
        elseif ismember(region_level, [4, 5])
            n_net = numel(net);
            sig_rt_mat = zeros(n_net, n_net);
            sig_rt_cor_mat = zeros(n_net, n_net);
            sig_lt_mat = zeros(n_net, n_net);
            sig_lt_cor_mat = zeros(n_net, n_net);
            ii = 1;
            for ni = 1:n_net
                sig_rt_mat(ni, ni:n_net) = h_rt(ii:(ii+n_net-ni));
                sig_rt_cor_mat(ni, ni:n_net) = h_rt_cor(ii:(ii+n_net-ni));
                sig_lt_mat(ni, ni:n_net) = h_lt(ii:(ii+n_net-ni));
                sig_lt_cor_mat(ni, ni:n_net) = h_lt_cor(ii:(ii+n_net-ni));
                ii = ii + n_net - ni + 1;
            end
            figure;
            subplot(1,2,1); imagesc(sig_rt_mat); xticks(1:n_net); xticklabels(net(1:n_net)); yticklabels(net);
            subplot(1,2,2); imagesc(sig_lt_mat); xticks(1:n_net); xticklabels(net(1:n_net)); yticklabels(net);
            
            Sig_rt_mat = cat(3, Sig_rt_mat, sig_rt_mat);
            Sig_rt_cor_mat = cat(3, Sig_rt_cor_mat, sig_rt_cor_mat);
            Sig_lt_mat = cat(3, Sig_lt_mat, sig_lt_mat);
            Sig_lt_cor_mat = cat(3, Sig_lt_cor_mat, sig_lt_cor_mat);
        end
        iii = iii + 1;
    end
    H_whole{compar,1} = Sig_rt_mat; H_whole_cor{compar,1} = Sig_rt_cor_mat;
    H_whole{compar,2} = Sig_lt_mat; H_whole_cor{compar,2} = Sig_lt_cor_mat;
end
%%
if region_level == 5
    H_whole_dmn = H_whole;
    H_whole_cor_dmn = H_whole_cor;
    save([root_path, 'results/graph_met/bct_metric_DS', str, '.mat'], 'H_whole_dmn', 'H_whole_cor_dmn', '-append')
end
%% save figure for edge level
if region_level == 2
    for compar = 1:length(x1_idx_cell)
        for ii = 1:4
            fig_name = ['set_',num2str(compar),'f_save', num2str(ii)];
            savefig(eval(fig_name), [root_path, 'results/Figures/DS',sprintf('%02d', mi), '_edge_ttest_',ft_name{x1_idx_cell{compar}(1)},'_',ft_name{x1_idx_cell{compar}(2)},'_fig', num2str(ii), '.fig'])
        end
    end
end