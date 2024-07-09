%% Leading Eigenvector_GhostAttractor
clear all
clc

% add the path of LEiDA toolbox
str = split(pwd, '/');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end

addpath(genpath([root_path, 'tools/LEiDA']))

%% bring timeseries from specific ROI including habenular, VTA, SNc
data_path = [root_path, 'data/RL.fanaticor.REML_wherr/'];
cd([data_path, '/Midbrain/SNc/'])

S_snc=dir('*.1D');
N_snc = {S_snc.name};
Nsubj = length(N_snc);

snc = [];   % result size = [108(=54*2) x 240]
for i = 1:Nsubj
	if ismember(i, [43, 44])
		continue
	else
		snc = [snc; load(N_snc{i})'];  % [1 x 240] L -> R
	end
end

cd([data_path, '/Midbrain/VTA/'])

S_vta=dir('*.1D');
N_vta = {S_vta.name};
Nsubj = length(N_vta);

vta = [];
for i = 1:Nsubj
	if ismember(i, [43, 44])
		continue
	else
		vta = [vta; load(N_vta{i})'];  % [1 x 240] L -> R
	end
end

cd([data_path, '/Midbrain/Habe/'])

S_habe=dir('*.1D');
N_habe = {S_habe.name};
Nsubj = length(N_habe);

habe = [];
for i = 1:Nsubj
	habe = [habe; load(N_habe{i})'];  % [1 x 240] L -> R
end

%%  Load data & calculate (first) leading eigenvector
cd([data_path, '/schaefer100-yeo7/MelbournSubcortex-S1/'])
S=dir('*.csv');
N = {S.name};
BOLD_tmp=load(N{2});
[~, Tmax]=size(BOLD_tmp);
Tmax = Tmax - 9;
n_Task = 1;
TR=1.5;
Nsubjs = length(N);
tian_roi = [1:5, 7:13, 15:16] + 100;
tian_roi = [1:100, tian_roi];
N_areas = length(tian_roi);

cd([data_path, '/schaefer100-yeo7/CIT168/'])

S_add1=dir('*.csv');
N_add1 = {S_add1.name};
add1_roi = [5:16, 21:32] + 100; 

cd([root_path, 'data/ANN/'])
S_add2=dir('*.mat');
N_add2 = {S_add2.name};
add2_roi = [1, 3, 6, 8, 10, 12, 14, 17, 20];
%%
N_seeds = N_areas+length(add1_roi)+length(add2_roi);  % 147
V1_all = zeros(Tmax*Nsubjs,N_seeds); % All leading eigenvectors
t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*(Tmax-2))
Time_all= zeros(Tmax*Nsubjs,1); % Vector that links each frame to a subject
FCD_eig=cell(Nsubjs, n_Task);
FCD_iFC=cell(Nsubjs, n_Task);
iFC_values=zeros(Tmax,(N_seeds*(N_seeds-1)/2));

ii=0;
%ts_all = [];
for si=1:Nsubjs
    if ismember(si, [22, 36, 44, 57])
        continue
    else
        disp(['Do subject : ', num2str(ii+1)])
        % Get the BOLD signals from subject si
        cd([data_path, '/schaefer100-yeo7/MelbournSubcortex-S1/'])
        BOLD = load(N{si});
        BOLD = BOLD(tian_roi,:);
        
        cd([data_path, '/schaefer100-yeo7/CIT168/'])
        BOLD_add1 = load(N_add1{si});
        BOLD_add1 = BOLD_add1(add1_roi,:);  % [24 x 240]
        
        BOLD_add1(4,:) = snc(ii * 2 + 2, :);    % Right
        BOLD_add1(9,:) = vta(ii * 2 + 2, :);   % Right
        BOLD_add1(11,:) = habe(ii * 2 + 2, :);  % Right
        BOLD_add1(16,:) = snc(ii * 2 + 1, :);    % Left
        BOLD_add1(21,:) = vta(ii * 2 + 1, :);    % Left
        BOLD_add1(23,:) = habe(ii * 2 + 1, :);   % Left
        BOLD = [BOLD; BOLD_add1];          % [138 x 240]
        
        cd([root_path, 'data/ANN/'])
        load(N_add2{si});
        BOLD_add2 = ts_AAN(add2_roi,:);  % [9 x 240]
        BOLD = [BOLD; BOLD_add2];          % [147 x 240]
        %ts_all = [ts_all, BOLD(:,6:Tmax+5)];
        
        ii = ii + 1;
        % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
        Phase_BOLD=zeros(N_seeds,Tmax+9);  % [147 x 240]
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_seeds
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
            %signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));  %filtered
            %Phase_BOLD(seed,:) = angle(hilbert(signal_filt)); %filtered
            Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:)));
        end
        
        %Calculate the Kuramoto Order Parameter
        OP=abs(sum(exp(1i*Phase_BOLD))/N_seeds);
        Metasta(si)=std(OP);
        Synchro(si)=mean(OP);
        
        %
        [V1previous,~]=eigs(cos(Phase_BOLD(:,1)-Phase_BOLD(:,1)'),1);   % [Nnodes x Nnodes] -> [Nnodes x 1]
        
        % Slide over time discarding the first and last epochs
        for t=6:Tmax+5  % 6:236
            
            %Calculate the Instantaneous FC (BOLD Phase Synchrony(Coherence))
            iFC=zeros(N_seeds);
            for n=1:N_seeds
                for p=1:N_seeds
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t)');
                end
            end
            
            % Get the leading eigenvector of Instantaneous BOLD Phases
            [V1,~]=eigs(cos(Phase_BOLD(:,t)-Phase_BOLD(:,t)'),1);
            
            % OPTION 1: half-switching
            if sum(V1)>0
                V1=-V1;
            end
            
            % OPTION 2: continous-switching
            distV1vsV1previous = 1-pdist([V1,V1previous]','cosine');
            distmirrorV1vsV1previous = 1-pdist([-V1,V1previous]','cosine');
            if distV1vsV1previous > distmirrorV1vsV1previous
                V1previous = V1;
            elseif distmirrorV1vsV1previous > distV1vsV1previous
                V1previous = -V1;
            end
            
            % Save V1 from all frames in all fMRI sessions
            t_all = t_all+1; % Update time
            Time_all(t_all) = si; % Information that at t_all, V1 corresponds to subject s in a given task
            iFC_values(t_all,:)=iFC(triu(ones(N_seeds),1)>0);
            % OPTION 1: half-switching
            Leading_Eig(t_all,:)=V1; %V1_all(t_all,:) = V1;
            
            % OPTION 2: continous-switching
            Leading_Eig_cont(t_all,:)=V1previous; %V1cont_all(t_all,:) = V1previous;
            
        end
        
        % Calculate the FCD(Cosine similarity of eigenvectors over time)
        
        for t1=(t_all - Tmax + 1):t_all
            eig1=squeeze(Leading_Eig(t1,:));
            iFC1=iFC_values(t1,:);
            for t2=(t_all - Tmax + 1):t_all
                eig2=squeeze(Leading_Eig(t2,:));
                iFC2=iFC_values(t2,:);
                FCD_eig{si,1}(t1,t2)=dot(eig1,eig2)/norm(eig1)/norm(eig2);
                FCD_iFC{si,1}(t1,t2)=dot(iFC1,iFC2)/norm(iFC1)/norm(iFC2);
            end
        end
    end
end
%%
n_net = 7;
load([root_path 'data/shaefer100_subcortex_',num2str(n_net),'system.mat']);

if n_net == 7
    net = {'VIS','SOM','DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SC'};
elseif n_net == 17
    net = {'VISCent', 'VISPer', 'SOMa','SOMb', 'DATa', 'DATb', 'VATa', 'VATb', 'LIMa', 'LIMb', 'ContA', 'ContB', 'ContC', 'DMNa', 'DMNb', 'DMNc', 'TP', 'SC'};
end

load([root_path, 'data/connectivity_sign.mat'], 'H')

Nnodes = size(Leading_Eig, 2);
labs = ones(Nnodes,1) * (n_net+1);
labs(1:length(lab)) = lab;
Leida_con = zeros(Nnodes, Nnodes, size(Leading_Eig, 1));

labs_dmn = ones(Nnodes,1) * (n_net+3);
labs_dmn(1:length(lab)) = lab;
labs_dmn(38:50) = [1 1 1 1 1 1 2 2 1 1 2 3 2] + 6;
labs_dmn(90:100) = [2 2 1 1 1 1 2 1 2 3 2] + 6;

inter = []; inter_dmn = [];
for ii = 1:size(Leading_Eig, 1)
    Leading_Eig_tmp = squeeze(Leading_Eig(ii, :));
    Leida_con_tmp = Leading_Eig_tmp' * Leading_Eig_tmp;
    Leida_con_tmp = rescale(Leida_con_tmp .* H);
    [inter_tmp, intra_tmp] = integration_recruitment(Leida_con_tmp, labs);
    [inter_dmn_tmp, intra_dmn_tmp] = integration_recruitment(Leida_con_tmp, labs_dmn);
    inter_tmp = inter_tmp - diag(diag(inter_tmp));
    inter_tmp = inter_tmp + diag(diag(intra_tmp));
    inter = cat(3, inter, inter_tmp(1:n_net, 1:n_net, :));  % (7, 7, tp*Nsubj)
    
    inter_dmn_tmp = inter_dmn_tmp - diag(diag(inter_dmn_tmp));
    inter_dmn_tmp = inter_dmn_tmp + diag(diag(intra_dmn_tmp));
    inter_dmn = cat(3, inter_dmn, inter_dmn_tmp(1:n_net + 2, 1:n_net + 2, :));
end

inter_vec = []; 
for ni = 1:n_net - 1
    inter_vec = [inter_vec; squeeze(inter(ni, ni:end,:))];
end
inter_vec = [inter_vec', squeeze(inter(n_net, n_net,:))];  % (28 x tp*Nsubj)

inter_dmn_vec = [];
for ni = 1:n_net + 2 - 1
    inter_dmn_vec = [inter_dmn_vec; squeeze(inter_dmn(ni, ni:end,:))];
end
inter_dmn_vec = [inter_dmn_vec', squeeze(inter_dmn(n_net + 2, n_net + 2,:))];  % (28 x tp*Nsubj)

%% save
check = exist([root_path, 'results/leida/'], 'dir');
if check == 0
    mkdir([root_path, 'results/leida/'])
end
save([root_path,'results/leida/Leida_eig.mat'], 'Leading_Eig', 'inter_vec') %, 'Leading_Eig_cont', 'Time_all', 'FCD_iFC', 'FCD_eig', '-v7.3');
save([root_path,'results/leida/Leida_eig.mat'], 'inter_dmn_vec', '-append')
