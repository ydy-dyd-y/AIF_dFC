clear
close all
clc

str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

n_net = 7;
thres = 0.01;
load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);

addpath(genpath([root_path,'tools/']));
load([root_path,'tools/fcn/surfinfo.mat'])
aa = niftiread([root_path,'tools/fcn/Schaefer2018_100Parcels_',num2str(n_net),'Networks_order.dlabel.nii']);
load([root_path, 'data/cm.mat'], 'cm')

aa = niftiread([root_path,'tools/fcn/Schaefer2018_100Parcels_',num2str(n_net),'Networks_order.dlabel.nii']);
aa = squeeze(aa);
ar = aa(32493:end); al = aa(1:32492);
cr = zeros(size(ar)); cl = zeros(size(al));
cmap = fcn_cmaphot;

if n_net == 7
    net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SC'};
elseif n_net == 17
    net = {'VISCent', 'VISPer', 'SOMa','SOMb', 'DATa', 'DATb', 'VATa', 'VATb', 'LIMa', 'LIMb', 'ContA', 'ContB', 'ContC', 'DMNa', 'DMNb', 'DMNc', 'TP', 'SC'};
end

Nnodes = length(lab_lr);
labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
lab = labs(1:100);
[gx,gy,idx] = grid_communities(lab); % BCT function

%%
try load([root_path, 'results/behav/hgf.mat'])
catch est_hgf_lore = cell(Nsubj, 2);
end

rm_sub1 = [20, 32];
rm_sub2 = [22, 36, 44, 57];
est_hgf_lore(rm_sub1,:) = [];
est_hgf_lore(rm_sub2,:) = [];

%% Loss task
start_i = 14;
tmp_si = [];
for si = 1:54
    tmp = [est_hgf_lore{si, 1}.u(start_i) == est_hgf_lore{si, 1}.y(start_i), est_hgf_lore{si, 1}.u(start_i+1) == est_hgf_lore{si, 1}.y(start_i+1), est_hgf_lore{si, 1}.u(start_i+2) == est_hgf_lore{si, 1}.y(start_i+2)];
    tmp_si = [tmp_si; tmp];
end

loss1_si1 = []; loss1_si2 = []; loss1_si3 = []; loss1_si4 = [];
for si = 1:54
    if sum(tmp_si(si,:) == [0 0 0]) == 3
        loss1_si1 = [loss1_si1, si];
    elseif sum(tmp_si(si,:) == [0 1 0]) == 3
        loss1_si2 = [loss1_si2, si];
    elseif sum(tmp_si(si,:) == [0 0 1]) == 3
        loss1_si3 = [loss1_si3, si];
    elseif sum(tmp_si(si,:) == [0 1 1]) == 3
        loss1_si4 = [loss1_si4, si];
    end
end

start_i = 21;
tmp_si = [];
for si = 1:54
    tmp = [est_hgf_lore{si, 1}.u(start_i) == est_hgf_lore{si, 1}.y(start_i), est_hgf_lore{si, 1}.u(start_i+1) == est_hgf_lore{si, 1}.y(start_i+1), est_hgf_lore{si, 1}.u(start_i+2) == est_hgf_lore{si, 1}.y(start_i+2), est_hgf_lore{si, 1}.u(start_i+3) == est_hgf_lore{si, 1}.y(start_i+3)];
    tmp_si = [tmp_si; tmp];
end

loss2_si1 = []; loss2_si2 = []; loss2_si3 = []; loss2_si4 = []; loss2_si5 = []; loss2_si6 = []; loss2_si7 = []; loss2_si8 = [];
for si = 1:54
    if sum(tmp_si(si,:) == [0 0 0 0]) == 4
        loss2_si1 = [loss2_si1, si];
    elseif sum(tmp_si(si,:) == [0 0 0 1]) == 4
        loss2_si2 = [loss2_si2, si];
    elseif sum(tmp_si(si,:) == [0 0 1 0]) == 4
        loss2_si3 = [loss2_si3, si];
    elseif sum(tmp_si(si,:) == [0 0 1 1]) == 4
        loss2_si4 = [loss2_si4, si];
    elseif sum(tmp_si(si,:) == [0 1 0 0]) == 4
        loss2_si5 = [loss2_si5, si];
    elseif sum(tmp_si(si,:) == [0 1 0 1]) == 4
        loss2_si6 = [loss2_si6, si];
    elseif sum(tmp_si(si,:) == [0 1 1 0]) == 4
        loss2_si7 = [loss2_si7, si];
    elseif sum(tmp_si(si,:) == [0 1 1 1]) == 4
        loss2_si8 = [loss2_si8, si];
    end
end

%% reward task
start_i = 21;
tmp_si = [];
for si = 1:54
    tmp = [est_hgf_lore{si, 2}.u(start_i) == est_hgf_lore{si, 2}.y(start_i), est_hgf_lore{si, 2}.u(start_i+1) == est_hgf_lore{si, 2}.y(start_i+1), est_hgf_lore{si, 2}.u(start_i+2) == est_hgf_lore{si, 2}.y(start_i+2)];
    tmp_si = [tmp_si; tmp];
end

reward1_si1 = []; reward1_si2 = []; reward1_si3 = []; reward1_si4 = [];
for si = 1:54
    if sum(tmp_si(si,:) == [0 0 0]) == 3
        reward1_si1 = [reward1_si1, si];
    elseif sum(tmp_si(si,:) == [0 1 0]) == 3
        reward1_si2 = [reward1_si2, si];
    elseif sum(tmp_si(si,:) == [0 0 1]) == 3
        reward1_si3 = [reward1_si3, si];
    elseif sum(tmp_si(si,:) == [0 1 1]) == 3
        reward1_si4 = [reward1_si4, si];
    end
end

start_i = 35;
tmp_si = [];
for si = 1:54
    tmp = [est_hgf_lore{si, 2}.u(start_i) == est_hgf_lore{si, 2}.y(start_i), est_hgf_lore{si, 2}.u(start_i+1) == est_hgf_lore{si, 2}.y(start_i+1), ...
        est_hgf_lore{si, 2}.u(start_i+2) == est_hgf_lore{si, 2}.y(start_i+2), est_hgf_lore{si, 2}.u(start_i+3) == est_hgf_lore{si, 2}.y(start_i+3), ...
        est_hgf_lore{si, 2}.u(start_i+4) == est_hgf_lore{si, 2}.y(start_i+4)];
    tmp_si = [tmp_si; tmp];
end

reward2_si1 = []; reward2_si2 = []; reward2_si3 = []; reward2_si4 = []; reward2_si5 = []; reward2_si6 = []; reward2_si7 = []; reward2_si8 = [];
reward2_si9 = []; reward2_si10 = []; reward2_si11 = []; reward2_si12 = []; reward2_si13 = []; reward2_si14 = []; reward2_si15 = []; reward2_si16 = [];
for si = 1:54
    if sum(tmp_si(si,:) == [0 0 0 0 0]) == 5
        reward2_si1 = [reward2_si1, si];
    elseif sum(tmp_si(si,:) == [0 1 0 0 0]) == 5
        reward2_si2 = [reward2_si2, si];
    elseif sum(tmp_si(si,:) == [0 0 1 0 0]) == 5
        reward2_si3 = [reward2_si3, si];
    elseif sum(tmp_si(si,:) == [0 1 1 0 0]) == 5
        reward2_si4 = [reward2_si4, si];
    elseif sum(tmp_si(si,:) == [0 0 0 1 0]) == 5
        reward2_si5 = [reward2_si5, si];
    elseif sum(tmp_si(si,:) == [0 1 0 1 0]) == 5
        reward2_si6 = [reward2_si6, si];
    elseif sum(tmp_si(si,:) == [0 0 1 1 0]) == 5
        reward2_si7 = [reward2_si7, si];
    elseif sum(tmp_si(si,:) == [0 1 1 1 0]) == 5
        reward2_si8 = [reward2_si8, si];

    elseif sum(tmp_si(si,:) == [0 0 0 0 1]) == 5
        reward2_si9 = [reward2_si9, si];
    elseif sum(tmp_si(si,:) == [0 1 0 0 1]) == 5
        reward2_si10 = [reward2_si10, si];
    elseif sum(tmp_si(si,:) == [0 0 1 0 1]) == 5
        reward2_si11 = [reward2_si11, si];
    elseif sum(tmp_si(si,:) == [0 1 1 0 1]) == 5
        reward2_si12 = [reward2_si12, si];
    elseif sum(tmp_si(si,:) == [0 0 0 1 1]) == 5
        reward2_si13 = [reward2_si13, si];
    elseif sum(tmp_si(si,:) == [0 1 0 1 1]) == 5
        reward2_si14 = [reward2_si14, si];
    elseif sum(tmp_si(si,:) == [0 0 1 1 1]) == 5
        reward2_si15 = [reward2_si15, si];
    elseif sum(tmp_si(si,:) == [0 1 1 1 1]) == 5
        reward2_si16 = [reward2_si16, si];
    end
end

confused_agent = union(union(union(union(union(loss1_si1, loss2_si5), reward1_si1), reward2_si1), reward2_si9), reward2_si3);

%% 
mi = 9;   % 7~9
str = sprintf('%02d', mi);
region_level = 4; behav_data_type = 1;
region_level_label = {'node', 'edge', 'state', 'system'}; method_label = {'act', 'dISC'};
data_type = [method_label{behav_data_type}, '_', region_level_label{region_level}, '_DS'];

tmp_name = [root_path, 'results/GLM/', data_type, num2str(str), '.mat'];
load(tmp_name, 'params_mass')

tmp_uc1 = squeeze(params_mass(:,:,5)); % uc2_1, (28 x 54)
tmp_uc2 = squeeze(params_mass(:,:,6)); % uc2_2
tmp_uc3 = squeeze(params_mass(:,:,7)); % uc2_3
tmp_uc4 = squeeze(params_mass(:,:,8)); % uc2_4
tmp_uc5 = squeeze(params_mass(:,:,9)); % uc2_5

A = cell({});
ii = 1;
for roi = 1:size(params_mass, 1)
    tmp1 = [tmp_uc1(roi,confused_agent)', tmp_uc2(roi,confused_agent)', tmp_uc3(roi,confused_agent)', tmp_uc4(roi,confused_agent)', tmp_uc5(roi,confused_agent)'];
    A{ii} = tmp1;
    tmp2 = zeros(length(confused_agent), 5);
    for si = 1:length(confused_agent)
        if ismember(confused_agent(si), loss1_si1)
            tmp2(si, 1) = 1;
        end
        if ismember(confused_agent(si), loss2_si5)
            tmp2(si, 2) = 1;
        end
        if ismember(confused_agent(si), reward1_si1)
            tmp2(si, 3) = 1;
        end
        if ismember(confused_agent(si), [reward2_si1, reward2_si9])
            tmp2(si, 4) = 1;
        end
        if ismember(confused_agent(si), reward2_si3)
            tmp2(si, 5) = 1;
        end
    end
    B = tmp2;
    ii = ii + 1;
end

save([root_path, 'data/',region_level_label{region_level},'_AB_DS',str,'.mat'], 'A', 'B')

%{
% checkerboard in Fig. 6c.
figure;
imagesc(B')
colormap([1 1 1; 0 0 0])
yticklabels([])
xticks(1:26)
xlabel('Subjects')
set(gcf, 'color', 'w')
set(gcf, 'Position', [680 237 647 741])
%}

uc_idx = 4;
eval(['tmp_uc = tmp_uc', num2str(uc_idx), ';'])
if uc_idx == 1
    uc_data = [tmp_uc(:,loss1_si1)'; tmp_uc(:,loss1_si3)'];
    grouplabels = [ones(length(loss1_si1), 1); zeros(length(loss1_si3), 1)];
elseif uc_idx == 2
    uc_data = [tmp_uc(:,loss2_si5)'; tmp_uc(:,loss2_si6)'];
    grouplabels = [ones(length(loss2_si5), 1); zeros(length(loss2_si6), 1)];
elseif uc_idx == 3
    uc_data = [tmp_uc(:,reward1_si1)'; tmp_uc(:,reward1_si3)'];
    grouplabels = [ones(length(reward1_si1), 1); zeros(length(reward1_si3), 1)];
elseif uc_idx == 4
    uc_data = [tmp_uc(:,[reward2_si1, reward2_si9])'; tmp_uc(:,[reward2_si3, reward2_si11])'];
    grouplabels = [ones(length([reward2_si1, reward2_si9]), 1); zeros(length([reward2_si3, reward2_si11]), 1)];
end

figure;
if ismember(region_level, [1, 2])
    H_rt = []; P_rt = []; H_lt = []; P_lt = [];
    for ri = 1:size(params_mass, 1)
        [h, p, ci, stats] = ttest2(uc_data(find(grouplabels == 1), ri), uc_data(find(grouplabels == 0), ri), 'Tail','right', 'Alpha', thres);
        H_rt = [H_rt, h]; P_rt = [P_rt, p];
        [h, p, ci, stats] = ttest2(uc_data(find(grouplabels == 0), ri), uc_data(find(grouplabels == 1), ri), 'Tail','right', 'Alpha', thres);
        H_lt = [H_lt, h]; P_lt = [P_lt, p];
    end
    if region_level == 1
        id = 0;
        while id < 2
            if id == 0
                data_tmp = H_rt;
            else
                data_tmp = H_lt;
            end
            cr(ar ~= 0) = data_tmp(ar(ar ~= 0));
            cl(al ~= 0) = data_tmp(al(al ~= 0));

            subplot(2, 4, 1 + id*2);
            if isfield(sr,'data')
                th = trisurf(sl.data{2}.data+1,sl.data{1}.data(:,1),sl.data{1}.data(:,2),sl.data{1}.data(:,3),cl);
            else
                th = trisurf(sl.faces,sl.vertices(:,1),sl.vertices(:,2),sl.vertices(:,3),cl);
            end
            set(th,'edgecolor','none'); axis image; %set(gca,'clim',cm_lim);
            view(gca,3);axis equal;axis off;view(-90,0);material dull;camlight headlight;lighting gouraud
            colormap(cm);%colorbar;
            subplot(2, 4, 6 + id*2);
            if isfield(sr,'data')
                th = trisurf(sr.data{2}.data+1,sr.data{1}.data(:,1),sr.data{1}.data(:,2),sr.data{1}.data(:,3),cr);
            else
                th = trisurf(sr.faces,sr.vertices(:,1),sr.vertices(:,2),sr.vertices(:,3),cr);
            end
            set(th,'edgecolor','none'); axis image; %set(gca,'clim',cm_lim);
            view(gca,3);axis equal;axis off;view(270,0);material dull;camlight headlight;lighting gouraud
            colormap(cm);

            subplot(2, 4, 2 + id*2);
            if isfield(sr,'data')
                th = trisurf(sr.data{2}.data+1,sr.data{1}.data(:,1),sr.data{1}.data(:,2),sr.data{1}.data(:,3),cr);
            else
                th = trisurf(sr.faces,sr.vertices(:,1),sr.vertices(:,2),sr.vertices(:,3),cr);
            end
            set(th,'edgecolor','none'); axis image; %set(gca,'clim',cm_lim);
            view(gca,3);axis equal;axis off;view(90,0);material dull;camlight headlight;lighting gouraud
            colormap(cm);
            subplot(2, 4, 5 + id*2);
            if isfield(sr,'data')
                th = trisurf(sl.data{2}.data+1,sl.data{1}.data(:,1),sl.data{1}.data(:,2),sl.data{1}.data(:,3),cl);
            else
                th = trisurf(sl.faces,sl.vertices(:,1),sl.vertices(:,2),sl.vertices(:,3),cl);
            end
            set(th,'edgecolor','none'); axis image; %set(gca,'clim',cm_lim);
            view(gca,3);axis equal;axis off;view(-270,0);material dull;camlight headlight;lighting gouraud
            colormap(cm);  
            set(gcf, 'color', 'w')
            id = id + 1;
        end
        sgtitle(['uc idx ', num2str(uc_idx)])
        f = gcf;
        f.Position = [622 199 1239 699];
    elseif region_level == 2
        ir = [];
        ic = [];
        for i = 1:Nnodes
            ir = [ir;repelem(i, Nnodes - i)'];
            ic = [ic;linspace(i+1, Nnodes, Nnodes - i)'];
        end
        beta_tmp_rt = zeros(Nnodes); beta_tmp_lt = zeros(Nnodes);
        for ii = 1:size(params_mass, 1)
            beta_tmp_rt(ir(ii), ic(ii)) = H_rt(ii);
            beta_tmp_rt(ic(ii), ir(ii)) = H_rt(ii);
            beta_tmp_lt(ir(ii), ic(ii)) = H_lt(ii);
            beta_tmp_lt(ic(ii), ir(ii)) = H_lt(ii);
        end
        subplot(1,2,1); imagesc(beta_tmp_rt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5); colorbar('southoutside'); title('wrong > right')
        subplot(1,2,2); imagesc(beta_tmp_lt(idx,idx)); axis square; hold on; plot(gx,gy,'k','linewidth',1.5); colorbar('southoutside'); title('right > wrong')
        sgtitle(['uc idx = ', num2str(uc_idx)])
    end
else
    for ri = 1:size(params_mass, 1)
        if region_level == 4
            if ri < 8
                subplot(7, 7, ri)
            elseif ri < 14
                subplot(7, 7, ri + 1)
            elseif ri < 19
                subplot(7, 7, ri + 3)
            elseif ri < 23
                subplot(7, 7, ri + 6)
            elseif ri < 26
                subplot(7, 7, ri + 10)
            elseif ri < 28
                subplot(7, 7, ri + 15)
            else
                subplot(7, 7, ri + 21)
            end
        end
        boxplot(uc_data(:,ri), grouplabels)
        [h, p, ci, stats] = ttest2(uc_data(find(grouplabels == 1), ri), uc_data(find(grouplabels == 0), ri), 'Tail','right');
        if h == 1
            fprintf('integration : %i\n', ri)
            fprintf('p-value : %.4f\n', p)
            title(sprintf('p-value : %.4f\n', p), 'Color','r')
        end
        if ismember(ri, [19, 21, 26])
            fprintf('integration : %i\n', ri)
            fprintf('p-value : %.4f\n', p)
        end
        [h, p, ci, stats] = ttest2(uc_data(find(grouplabels == 1), ri), uc_data(find(grouplabels == 0), ri), 'Tail','left');
        if h == 1
            fprintf('integration : %i\n', ri)
            fprintf('p-value : %.4f\n', p)
            title(sprintf('p-value : %.4f\n', p), 'Color','b')
        end
        if ismember(ri, [19, 21, 26])
            fprintf('integration : %i\n', ri)
            fprintf('p-value : %.4f\n', p)
        end

    end
    sgtitle(['Uncertain feedback2 ', num2str(uc_idx)])
end