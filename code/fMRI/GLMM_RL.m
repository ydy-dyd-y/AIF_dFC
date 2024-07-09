clear
clc
close all
% Graph Laplacian Mixture Model
str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end
addpath(genpath([root_path, 'tools/']))

%%
Nsubj = 54;
load([root_path,'data/ts_all.mat'], ['ts_all_',num2str(Nsubj),'sub'])
eval(['ts_all = ts_all_', num2str(Nsubj), 'sub;';])
tp = size(ts_all, 1)/Nsubj;

% Defining hyperparameters
k = 4;
delta = 0.3;
iterations = 200;
iter = 6;
spread = 0.1;
regul = 0.15;
ncores = 8;

Ls_iter = cell(iter,1);
mus_iter = zeros(size(ts_all,2), k, iter);
gamma_hats_iter = zeros(size(ts_all,1), k, iter);
W_iter = cell(iter,1);

idx_id = [];
idx_task = [];
for si = 1:Nsubj
    idx_id = [idx_id; si*ones(tp, 1)];
    idx_task = [idx_task; zeros(7,1); ones(107,1)*(-1); zeros(6,1); ones(107,1); zeros(4,1)];
end

try
    delete(gcp('nocreate'));
    pool = parpool('local',ncores);

    parfor it = 1:iter
        [Ls, gamma_hats, mus, log_likelihood, W] = glmm_matlab(ts_all, iterations, k, sign, spread,regul, delta);
        Ls_iter{it} = Ls;
        mus_iter(:,:,it) = mus;
        gamma_hats_iter(:,:,it) = gamma_hats;
        W_iter{it} = W;
    end
catch
    for it = 1:iter
        [Ls, gamma_hats, mus, log_likelihood, W] = glmm_matlab(ts_all, iterations, k, sign, spread,regul, delta);
        Ls_iter{it} = Ls;
        mus_iter(:,:,it) = mus;
        gamma_hats_iter(:,:,it) = gamma_hats;
        W_iter{it} = W;
    end
end

save([root_path, 'results/glmm/glmm_k' num2str(k) '_' num2str(delta) '.mat'], ...
    'gamma_hats_iter', 'mus_iter', 'Ls_iter', 'W_iter', 'ts_all', 'Nsubj', 'tp')

iter_tmp = 5;
gamma_hats = squeeze(gamma_hats_iter(:,:,iter_tmp));
[~, partitions] = max(gamma_hats, [], 2);  % [tp * Nsubj x 1]

% generate idx_all.mat 
% column 1 : subject id, column2 : task type(0:resting, -1:loss, 1:reward), column 3 : which state 
idx_all = [idx_id, idx_task, partitions];

var_name = ['idx_all_K',num2str(k),'_delta',num2str(delta * 10),'_iter',num2str(iter_tmp)];
eval([var_name ' = idx_all;']);

list = {dir([root_path, 'data']).name}; tf = contains(list, 'idx_all.mat'); 
if isempty(find(tf))
    save([root_path, 'data/idx_all.mat'], var_name)
else
    save([root_path, 'data/idx_all.mat'], var_name, '-append')  % already run "utility/indexing.m"
end

%% plotting the result states. 
% Fig. 3b~c
K = 4; delta = 0.3; iter_tmp = 5;
load([root_path, 'results/glmm/glmm_k',num2str(K),'_',num2str(delta),'.mat']);

W = W_iter{iter_tmp};
Ls = Ls_iter{iter_tmp};
Nnodes = size(W, 1);
A = zeros(Nnodes, Nnodes, K);
for k = 1: K 
    tmp = Ls(:,:,k);
    A(:,:,k) = diag(diag(tmp)) - tmp;  % this A has the same number of edges that 
                                       % have negative values with W estimated from sign = 1.
                                       % and the minimum value is also same.
end
mus = squeeze(mus_iter(:,:,iter_tmp));
gamma_hats = squeeze(gamma_hats_iter(:,:,iter_tmp));
gamma_hats_tmp = reshape(gamma_hats, [tp, Nsubj, K]);

% brain space plotting figures
n_net = 7;
load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);

% load 32k surfaces
load([root_path,'tools/fcn/surfinfo.mat'])
aa = niftiread([root_path,'tools/fcn/Schaefer2018_100Parcels_',num2str(n_net),'Networks_order.dlabel.nii']);
load([root_path, 'data/cm.mat'], 'cm')
cmap = fcn_cmaphot;
cmap_trim = cmap(65:end,:);

if n_net == 7
    net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SC'};
elseif n_net == 17
    net = {'VISCent', 'VISPer', 'SOMa','SOMb', 'DATa', 'DATb', 'VATa', 'VATb', 'LIMa', 'LIMb', 'ContA', 'ContB', 'ContC', 'DMNa', 'DMNb', 'DMNc', 'TP', 'SC'};
end

aa = squeeze(aa);
al = aa(1:32492);
ar = aa(32493:end);

labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
lab = labs;
lab(101:end) = [(2:8)'; (2:8)'; (10:2:32)'; (9:2:31)'; (33:41)'] + (n_net - 1);  % length of gx increase from 48 -> 210
[gx,gy,idx] = grid_communities(lab); % BCT function

subLabels = {'R-HIP','L-HIP','R-AMY','L-AMY','R-pTHAL','L-pTHAL','R-aTHAL','L-aTHAL', 'R-NAc','L-NAc','R-PUT','L-PUT','R-CAU','L-CAU',...
    'R-GPe','L-GPe','R-GPi','L-GPi','R-SNc','L-SNc','R-red','L-red','R-SNr','L-SNr','R-PBPN','L-PBPN',...
	'R-VTA','L-VTA', 'R-VP','L-VP','R-HABN','L-HABN','R-hypoTHAL','L-hypoTHAL','R-MM','L-MM', 'R-STN','L-STN', ...
    'DR', 'LC', 'LDTg', 'MnR', 'mRt', 'PAG', 'PBC', 'PnO', 'PTg'};

for k = 1:K
    mu = mus(:,k);

    eval(['f', num2str(k), '_Act', ' = figure';])
    subplot(4,2,[1,2]);
    boxplot(mu,labs,'labels',net,'labelorientation','horizontal'); 
    %title('Activation pattern per system(\mu)')
    
    mu_rank = mu(1:100); 
    cr = ones(size(ar)) * min(mu_rank);
    cl = ones(size(al)) * min(mu_rank);
    cr(ar ~= 0) = mu_rank(ar(ar ~= 0));
    cl(al ~= 0) = mu_rank(al(al ~= 0));
    %clim_val = [(-1)*max(abs(mu_rank)) max(abs(mu_rank))]; 
    subplot(4,2,3);
    if isfield(sr,'data') 
        th = trisurf(sl.data{2}.data+1,sl.data{1}.data(:,1),sl.data{1}.data(:,2),sl.data{1}.data(:,3),cl);
    else
        th = trisurf(sl.faces,sl.vertices(:,1),sl.vertices(:,2),sl.vertices(:,3),cl);
    end
    set(th,'edgecolor','none'); axis image; set(gca,'clim',[min(mu_rank),max(mu_rank)]);
    view(gca,3);axis equal;axis off;view(-90,0);material dull;camlight headlight;lighting gouraud
    
    subplot(4,2,4);
    if isfield(sr,'data') 
        th = trisurf(sr.data{2}.data+1,sr.data{1}.data(:,1),sr.data{1}.data(:,2),sr.data{1}.data(:,3),cr);
    else
        th = trisurf(sr.faces,sr.vertices(:,1),sr.vertices(:,2),sr.vertices(:,3),cr);
    end    
    set(th,'edgecolor','none'); axis image; set(gca,'clim',[min(mu_rank),max(mu_rank)]);
    view(gca,3);axis equal;axis off;view(90,0);material dull;camlight headlight;lighting gouraud
    
    subplot(4, 2, 5);
    if isfield(sr,'data') 
        th = trisurf(sl.data{2}.data+1,sl.data{1}.data(:,1),sl.data{1}.data(:,2),sl.data{1}.data(:,3),cl);
    else
        th = trisurf(sl.faces,sl.vertices(:,1),sl.vertices(:,2),sl.vertices(:,3),cl);
    end     
    set(th,'edgecolor','none'); axis image; set(gca,'clim',[min(mu_rank),max(mu_rank)]);
    view(gca,3);axis equal;axis off;view(-270,0);material dull;camlight headlight;lighting gouraud
    
    subplot(4, 2, 6);
    if isfield(sr,'data') 
        th = trisurf(sr.data{2}.data+1,sr.data{1}.data(:,1),sr.data{1}.data(:,2),sr.data{1}.data(:,3),cr);
    else
        th = trisurf(sr.faces,sr.vertices(:,1),sr.vertices(:,2),sr.vertices(:,3),cr);
    end
    set(th,'edgecolor','none'); axis image; set(gca,'clim',[min(mu_rank),max(mu_rank)]);
    view(gca,3);axis equal;axis off;view(270,0);material dull;camlight headlight;lighting gouraud
    
    subplot(4,2,[7,8])
    if isfield(sr,'data') 
        th = trisurf(sr.data{2}.data+1,sr.data{1}.data(:,1),sr.data{1}.data(:,2),sr.data{1}.data(:,3),cr);
    else
        th = trisurf(sr.faces,sr.vertices(:,1),sr.vertices(:,2),sr.vertices(:,3),cr);
    end
    set(th,'edgecolor','none'); axis image; set(gca,'clim',[min(mu_rank),max(mu_rank)]);
    view(gca,3);axis equal;axis off;view(270,0);material dull;camlight headlight;lighting gouraud
    %colormap(cmap);colorbar('southoutside');
    colormap(cmap);colorbar('southoutside'); 
    
    sgtitle(sprintf('State=%d ', k))
    set(gcf,'color','w', 'position', [50, 50, 800, 1500]);
    
    eval(['f', num2str(k), '_Adj', ' = figure';])
    %subplot(4, 2, [7, 8])
    W_tmp = A(:,:,k);
    %W_tmp = W(:,:,k);
    imagesc(W_tmp(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5); 

    if n_net == 7
        for i = 1:numel(net)-1
            labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 12);
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
    elseif n_net == 17
        for i = 1:numel(net)
            if mod(i,2) ==1
                labelhandle = text(mean([gx(6*i-1), gx(6*i-2)]), -2.5, net{i}, 'FontSize', 7.5);
                labelhandle.HorizontalAlignment = 'center';
            else
                labelhandle = text(-5.5, mean([gy(6*i-5), gy(6*i-4)]), net{i}, 'FontSize', 7.5);
                labelhandle.HorizontalAlignment = 'center';
            end
            yticks([])
            xticks([])
        end
    end
    colormap(cmap_trim);colorbar('southoutside');caxis([0 1])
    set(gcf,'color','w', 'position', [50, 50, 1800, 1800]);
end

%% polygonal plot of similarity with RSNs.
[~, ~, ~, tmp1] = NAME_CLUSTERS_ANGLE8(mus);
net8angle_pos = tmp1(:,1:n_net+1);
net8angle_neg = tmp1(:,n_net+2:(n_net+1)*2);
opt_axes.Labels = net;

opt_area.err        = 'std';
opt_area.FaceAlpha  = 0;
opt_area.Color      = [1 1 1; 1 1 1]./255;
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Labels    = false;
opt_lines.Legend    = {'neg','pos'};
opt_lines.Color     = [ 52 148 186; 236 112  22]./255;

f = figure;
set(f, 'position', [100, 100, 1700, 800])
for k = 1:K
    subplot(1,fix(K/1)+1,k);
    d_ex = [net8angle_neg(k,:); net8angle_pos(k,:)];
    d_ex = d_ex';
    d_ex = ( d_ex > 0 ).* d_ex;
    data = d_ex;
    polygonplot(data,opt_axes,opt_lines,opt_area);
end
%% save
savedir = [root_path, 'results/glmm/fig/'];
for ki = 1:K
    savefile = fullfile(savedir,['Fig.3_State',num2str(ki),'.pdf']);
    eval(['ax = f',num2str(ki),'_Act'])
    exportgraphics(ax, savefile) 
    eval(['ax = f',num2str(ki),'_Adj'])
    exportgraphics(ax, savefile,'Append',true) % 'Append' works after R2021b
end

%% plotting Fig. 3a
var_name = ['idx_all_K',num2str(k),'_delta',num2str(delta * 10),'_iter',num2str(iter_tmp)];
load([root_path, 'data/idx_all.mat'], var_name)
eval(['idx_all = ' var_name ';']);
ts_state = zeros(Nnodes, Nsubj, K);
for si = 1:Nsubj
    for ki = 1:K
        tmp_ts = [];
        for tt = 1:size(ts_all,1)
            if (idx_all(tt, 3) == ki) && (idx_all(tt, 1) == si)
                tmp_ts = [tmp_ts; ts_all(tt, :)];
                ts_state(:, si, ki) = mean(tmp_ts);
            end
        end
    end
end

load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);
labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
lab = labs;
lab(101:end) = [(2:8)'; (2:8)'; (10:2:32)'; (9:2:31)'; (33:41)'] + (n_net - 1);   % length of gx increase from 48 -> 210
[gx,gy,idx] = grid_communities(lab); % BCT function

subLabels = {'R-HIP','L-HIP','R-AMY','L-AMY','R-pTHAL','L-pTHAL','R-aTHAL','L-aTHAL', 'R-NAc','L-NAc','R-PUT','L-PUT','R-CAU','L-CAU',...
    'R-GPe','L-GPe','R-GPi','L-GPi','R-SNc','L-SNc','R-red','L-red','R-SNr','L-SNr','R-PBPN','L-PBPN',...
	'R-VTA','L-VTA', 'R-VP','L-VP','R-HABN','L-HABN','R-hypoTHAL','L-hypoTHAL','R-MM','L-MM', 'R-STN','L-STN', ...
    'DR', 'LC', 'LDTg', 'MnR', 'mRt', 'PAG', 'PBC', 'PnO', 'PTg'};

mus = mus(idx, :); % (147 x 4)
ts_state = ts_state(idx, :, :);
ts_sub_state = ts_state(101:end, :, :);

%% ANOVA
P = [];
thres = 0.05;
lt = [];
C = []; pair = cell({});
for roi = 1:size(ts_sub_state, 1)
    [p, tbl, stats] = anova1(squeeze(ts_sub_state(roi, :, :)), {'state1','state2','state3','state4'}, 'off');
    P = [P, p];
    if p < thres
        [c,m,h,gnames] = multcompare(stats);
        tmp_pair = find(c(:,end) < thres);
        if isempty(tmp_pair)
            tmp_pair = find(c(:,end) < 0.05);
        end
        lt = [lt, length(tmp_pair)];
        C = cat(3, C, c);
        pair{length(lt)} = tmp_pair;
    end
end

H = find(P < thres);
RED = cell({}); BLUE = cell({});
for ii = 1:length(H)
    c = squeeze(C(:,:,ii));
    if lt(ii) == 1
        if c(pair{ii}, 4) > 0
            RED{ii} = c(pair{ii}, 1); BLUE{ii} = c(pair{ii}, 2);
        else
            RED{ii} = c(pair{ii}, 2); BLUE{ii} = c(pair{ii}, 1);
        end
    elseif lt(ii) == 2
        if c(pair{ii}(1) ,(c(pair{ii}(1), 4) < 0) + 1) == c(pair{ii}(2) ,(c(pair{ii}(2), 4) < 0) + 1)
            RED{ii} = c(pair{ii}(1) ,(c(pair{ii}(1), 4) < 0) + 1);
            [~, min_idx] = max([abs(c(pair{ii}(1), 4)), abs(c(pair{ii}(2), 4))]);
            BLUE{ii} = c(pair{ii}(min_idx), (c(pair{ii}(min_idx), 4) < 0));
        else
            [~, max_idx] = max([abs(c(pair{ii}(1), 4)), abs(c(pair{ii}(2), 4))]);
            RED{ii} = c(pair{ii}(max_idx), (c(pair{ii}(max_idx), 4) < 0) + 1);
            BLUE{ii} = c(pair{ii}(1) ,(c(pair{ii}(1), 4) > 0) + 1);
        end
    elseif lt(ii) == 3
        if c(pair{ii}(1), (c(pair{ii}(1), 4) < 0) + 1) == c(pair{ii}(2), (c(pair{ii}(2), 4) < 0) + 1) && c(pair{ii}(1) ,(c(pair{ii}(1), 4) < 0) + 1) == c(pair{ii}(3) ,(c(pair{ii}(3), 4) < 0) + 1)
            RED{ii} = c(pair{ii}(1) ,(c(pair{ii}(1), 4) < 0) + 1);
        elseif c(pair{ii}(1), (c(pair{ii}(1), 4) > 0) + 1) == c(pair{ii}(2), (c(pair{ii}(2), 4) > 0) + 1) && c(pair{ii}(1),(c(pair{ii}(1), 4) > 0) + 1) == c(pair{ii}(3) ,(c(pair{ii}(3), 4) > 0) + 1)
            BLUE{ii} = c(pair{ii}(1) ,(c(pair{ii}(1), 4) > 0) + 1);
        else
            [~, max_idx] = max(c(:,5));
            RED{ii} = c(max_idx, (c(max_idx, 4) < 0) + 1);
            [~, min_idx] = min(c(:,5));
            BLUE{ii} = c(min_idx, (c(min_idx, 4) < 0));
        end
    elseif lt(ii) == 4
        RED{ii} = unique([c(pair{ii}(1) ,(c(pair{ii}(1), 4) < 0) + 1), c(pair{ii}(2) ,(c(pair{ii}(2), 4) < 0) + 1), c(pair{ii}(3) ,(c(pair{ii}(3), 4) < 0) + 1), c(pair{ii}(4) ,(c(pair{ii}(4), 4) < 0) + 1)]);
        BLUE{ii} = unique([c(pair{ii}(1) ,(c(pair{ii}(1), 4) > 0) + 1), c(pair{ii}(2) ,(c(pair{ii}(2), 4) > 0) + 1), c(pair{ii}(3) ,(c(pair{ii}(3), 4) > 0) + 1), c(pair{ii}(4) ,(c(pair{ii}(4), 4) > 0) + 1)]);
    else
        [~, max_idx] = max(c(:,5));
        RED{ii} = c(max_idx, (c(max_idx, 4) < 0) + 1);
        [~, min_idx] = min(c(:,5));
        BLUE{ii} = c(min_idx, (c(min_idx, 4) < 0));
    end
end

for ik = 1:K
    red = []; blue = [];
    sub_mu = mus(101:138,ik);
    for ii = 1:numel(RED)
        if ismember(ik, RED{ii})
            red = [red, H(ii)];
        end
        if ismember(ik, BLUE{ii})
            blue = [blue, H(ii)];
        end
    end
    figure;
    for isc = 1:3
        subplot(3, 1, isc); 
        var_name = ['b', num2str(isc)];
        if isc == 1
            tmp = 1:14;
            y_lim = [-0.25, 0.25];
        elseif isc == 2
            tmp = 15:26;
            y_lim = [-0.05, 0.05];
        else
            tmp = 27:38;
            y_lim = [-0.1, 0.1];
        end
        red_tmp = find(ismember(tmp, red));
        blue_tmp = find(ismember(tmp, blue));
        eval([var_name, ' = bar(sub_mu(tmp), "k");'])
        eval([var_name, '.FaceColor = "flat";'])
        eval([var_name, '.CData(red_tmp,:) = repmat([0.8350 0.0780 0.1240], length(red_tmp), 1);'])
        eval([var_name, '.CData(blue_tmp,:) = repmat([0 0.4470 0.7410], length(blue_tmp), 1);'])
        xticks(1:length(tmp))
        xticklabels(subLabels(tmp));
        ylim(y_lim)
        t = sgtitle(['State', num2str(ik)]);
        t.FontSize = 14;        
    end
    set(gcf, 'color', 'w')
end
