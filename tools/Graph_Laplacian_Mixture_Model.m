%% Graph Laplacian Mixture Mode
addpath('/APP/NeuroScience/NeuroscienceToolbox/GLMM')
addpath(genpath('/APP/NeuroScience/NeuroscienceToolbox/gspbox'))
addpath(genpath('/APP/NeuroScience/NeuroscienceToolbox/unlocbox'))

%% sources
% https://github.com/Hermina/GLMM/blob/master/glmm_matlab.m
% https://github.com/epfl-lts2/unlocbox/releases/tag/1.7.5
% https://github.com/epfl-lts2/unlocbox
% https://github.com/epfl-lts2/gspbox
% normalized mutual information: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/1c38ca31-840b-4491-ac07-bd3ef7b18ad3/ccec48ae-2dc7-451d-a705-3790374f6b5f/previews/FSLib_v6.2.1_2018/lib/nmi.m/index.html

%% articles
% https://www.sciencedirect.com/science/article/pii/S1053811922001665
% (Dynamics of functional network organization through graph mixture learning)
% appendix: https://ars.els-cdn.com/content/image/1-s2.0-S1053811922001665-mmc1.pdf

% https://www.sciencedirect.com/science/article/pii/S1053811922001008?via%3Dihub#sec0024


%%
cd '/u5/Dense/ds002674/derivatives/Results/meants_corrmap/rfMRI_REST1_PA_fmap/timeseries/schaefer100-yeo7/MelbournSubcortex-S1/'

%ts=csvread('YAD10004.afni.schaefer100-yeo7_MelbournSubcortex-S1.csv');
%ts=csvread('ses-01.rfMRI_REST1_PA_fmap.schaefer100-yeo7_MelbournSubcortex-S1.csv');
S=dir('*.csv');
N = {S.name};
ts_all = [];
TR=0.72;
HighPass=0.008; %1/128 sec
LowPass=0.1;  % TR=0.72

for i=1:length(N)
%for i=1:6
  ts=csvread(N{i});
    
  mat_size=size(ts,1); % number of ROIs
  tp=size(ts,2);  % number of time points

% band pass filtering [0.01 0.1]
% https://m.blog.naver.com/PostView.naver?isHttpsRedirect=true&blogId=coolstu&logNo=130115108264
% https://github.com/DCAN-Labs/dcan_bold_processing/blob/master/matlab_code/filtered_movement_regressors.m
% https://github.com/lindenmp/rs-fMRI/blob/master/func/ButterFilt.m  (4th
% order Butter worth bandpass filtering
  ftype='bandpass';
  Fs=1/TR; % Fs = 1/TR(s)
  Fn=Fs/2; % Fn: Nyquist frequency = sampling freqency(Fs)/2 
% Bandpass cutoffs
  Wn=[HighPass LowPass];
% Build filter
  FiltOrder = 2;
  [bs,as] = butter(FiltOrder,Wn/Fn,ftype);  % filter order = 2 : ?4th order butterworth bandpass filtering
  ts_tp=ts';  % transpose
  %ts_bp=filter(bs,as,ts_tp);  % coarser than filtfilt
  ts_bp=filtfilt(bs,as,ts_tp);

% get rid of first & last time points (2/5 hz of band width  * TR)
% Chumin et al. (2022) Cortico-subcortical interactions in overlapping communities of edge functional connectivity
% https://www.sciencedirect.com/science/article/pii/S1053811922001008?via%3Dihub#bib0045
% https://github.com/brain-networks/edge-centric_demo

% Motion artifacts had been already performed during preprocessing steps
% Nilearn signal.clean, which removes confounds orthogonal to the temporal filters (Lindquist et al., 2019)

 half_bw=(Wn(2)-Wn(1))*2/5; 
 %half_bw=(Wn(2)-Wn(1))/2;
 half_bw_sec=1/half_bw;
 tp_remove=double(int8(half_bw_sec/TR));
 ts_bp_removed_tp=ts_bp(tp_remove:tp-tp_remove,:);

  % zscore time series
  %z = zscore(ts');
  z = zscore(ts_bp_removed_tp);

  ts_all = [ts_all; z];  % concatenate all subjects' ts

end


%% "train a glmm on data ts"
%% find optimal hyper-parameters (k, delta, gamma)
% The optimization of the number of clusters (𝐾) is also guided by the consistency of the
%  obtained clustering probabilities with the task paradigm. we evaluate the results by iteratively changing 𝐾 according
%  to the concordance of the estimated 𝜸 with respect to the experimental task paradigm
%  See Supplementary Material for further information (Appendix A.2.1, Fig. A2).
% On the other hand, RS data lacks a ground truth, so the number of
%  clusters 𝐾 has been chosen based on the optimized silhouette measure
%  and the consensus clustering procedure Monti et al. (2003), which is a
%  resampling-based method for optimal class discovery. 𝐾 has eventually
%  been varied according to the procedure mentioned above, in order to
%  capture multiple networks. In practice, changing the optimal number of
%  clusters does not seem to affect the final estimation, instead, it opens
%  the possibility of inferring more or fewer networks

% grid search (log_likelihood? Nop! cosine similarity)for 3 parameters: k, theta, delta
% in glmm_matlab: delta is 2 as default (5 ~ 15.5, 1.5 interval)
% k is needed to put as an input (2 ~ 11).
% in case of theta, norm_par (line 6: 0.65 ~ 2, 0.15 interval) is changed
% with function (theta = mean(Z(:))/norm_par) at line 53 of glmm_matlab.m

% in rfMR case, maximizing the silhouette is also applied. (suppl. A.12)
% silhouette(i) = b(i) - a(i) / max[a(i),b(i)]
% : where i is the single data point, a is the average
% distance between i and all other data points in the same cluster and b
% the smallest average distance of i to all points in any other cluster

% performance is cosine similarity (Figure A1 in suppl.)
% The optimization is implemented focusing solely on the activation
% patterns or mean activation (mu). On task fMRI we split the data in 2
% groups, run the GLMM algorithm and match the activation patterns of the
% two groups with the Hungarian algorithm. This is done in a grid-search
% varying the three hyper-parameters. The process has been bootstrapped to
% randomize the splits and get standard deviations of the estimates
% Eventually, we maximize the worst cosine similarity between the matched
% activation patterns (Figure A1: The worst similarity between two splits. the values of the metric with its standard deviations 
% in function of one of the three hyper-parameters, by fixing the other two
% to their optimal)
% 1) split into 2 groups, 2) run with certain value of a hyper-parameter in
% each group and get mus, 3) bootstrapping: run step 1-2 again again and
% get mus'distribution and its sd, 4) run step 1-3 with different value of
% a hyper-parameter, 5) run step 1-4 with other 2 hyper-parameters by
% fixing the remain two


%%

% k = 3 ; 
spread = 0.1;
regul = 0.1

% three parameters
% change value of a parameter by fixing the remain two
% add boostraping code to divide ts_all into two group randomly in following for-loop

%for k=2:9
k=4; %number of clustersv% need to optimise
norm_par = 1.5; %0.3; %1.2;
delta=2; % if want to change the value of delta, please make disable the line 2 in glmm_matlab.m

iterations = 200;
[Ls, gamma_hats, mus, log_likelihood] = glmm_matlab(ts_all, iterations,k,spread,regul,norm_par);
%disp('Training done')

% distance in activation pattern between two group for task fMR
% https://stackoverflow.com/questions/57187941/how-to-calculate-cosine-similarity-between-two-frequency-vectors-in-matlab
%Cos(x, y) = x . y / ||x|| * ||y||
%cosSim = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
%cosSim = (a(:).'*b(:))/sqrt(sum(a.^2)*sum(b.^2));
%dist=1-pdist([mus(3:4,:),mus(1:2,:)],'cosine');

% maximize the silhouette, defined as clustering metric for resting fMR
%clustev_cosine = evalclusters(Leading_Eig,'kmeans','silhouette','KList',2:19,'Distance','cosine');

%end


%%
disp(sum(gamma_hats,1)); % probability estimates
area(gamma_hats(1:745,:),'DisplayName','gamma_hats(1:745,:)') % the first subject's gamma hat ################# 231 -> 745
imagesc(Ls(:,:,1)); % 1st graph laplacian
imagesc(Ls(:,:,2)); % 2nd graph laplacian
% mus % activation patterns
plot(mus,'DisplayName','mus')


%%
%tp_remain=tp-tp_remove;
%gamma_hats_each=reshape(gamma_hats,[tp_remain,k,length(S)]);
%gamma_hats_each=reshape(gamma_hats',[231,4,54]);
%gamma_hats_each=reshape(gamma_hats,54,[],size(gamma_hats,2));
gamma_hats_each=reshape(gamma_hats,length(S),[],size(gamma_hats,2));
gamma_hats_mean=squeeze(mean(gamma_hats_each));
area(gamma_hats_mean)

%% "generate graphs"
n = size(ts_all,2); %15; %graph size

m = size(ts_all,1); %150; %number of signals

zero_thresh = 10e-4;

g(k) = gsp_erdos_renyi(n,0.7);

for i = 1:k
    while(1)
    	g(i) = gsp_erdos_renyi(n, 0.7);
        eigs = sort(eig(g(i).L));
        if (eigs(2) > zero_thresh) %ensuring graphs are connected
            break;
        end
    end
end

gamma = rand([m,1]);
gamma_cut = zeros(m,k);
dist = 0.5;
%p = [0, 0.2, 0.6, 1];
p = 0:1/k:1;
y = zeros(m,n);
true_y = zeros(m,n,k);
center = zeros(n,k);
gauss = zeros(n, n, k);
Lap = zeros(n, n, k);

for i=1:k
    gc = pinv(full(g(i).L));
    gauss(:,:,i) = (gc +gc')/2;
    Lap(:,:,i) = full(g(i).L);
    % generate centers, responsibilities and data
    center(:,i) = dist * randn([n,1]);
    center(:,i) = center(:,i) - mean(center(:,i));
    gamma_cut(p(i)<gamma & gamma<=p(i+1), i) = 1;
    true_y(:,:,i) = squeeze(gamma_cut(:,i)).*mvnrnd(center(:,i),gauss(:,:,i),m);
    y = y + true_y(:,:,i);
end

imagesc(Lap(:,:,2));
%%
% we need nmi function (normalized mutual information)
% https://kr.mathworks.com/matlabcentral/fileexchange/29047-normalized-mutual-information
% /APP/NeuroScience/NeuroscienceToolbox/nmi.m
addpath('/home/dayoung/program/')

[identify, precision, recall,  f, cl_errors] = identify_and_compare(Ls, Lap, gamma_hats, gamma_cut, k)

