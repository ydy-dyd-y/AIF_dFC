clear
clc
% Graph Laplacian Mixture Model
str = split(pwd, '\');
root_path = [];
for i = 1:numel(str)-2
    root_path = [root_path, str{i}, '/'];
end
addpath(genpath([root_path, 'tools/']))

%%
Nsubj = 54;
load([root_path,'data/ts_all.mat'], ['ts_all_',num2str(Nsubj),'sub'])
eval(['ts_all = ts_all_', num2str(Nsubj), 'sub';])
tp = 231;

%% Defining the range of hyperparameters
K = 3:8;
Delta = [0.3 0.6 0.9 1.2 1.5];
spread = 0.1;
regul = 0.15;
iterations = 200;
iter = 25;

inK = 1;
CD_MAX = zeros(iter,1);
cost = zeros(length(Delta),length(K));
statek_name = cell({});
delta_name = cell({});
ncores = 8;


for k = K
    statek_name{1,k} = num2str(k);
    inD =1;
    for delta = Delta
        delta_name{1,inD} = ['delta = ',num2str(delta)];
        % split the group data into 2 group 
        try
            delete(gcp('nocreate'));
            pool = parpool('local',ncores);

            parfor it = 1:iter
                idx = randperm(Nsubj);
                ts_1 = [];
                ts_2 = [];
                for a = idx(1:Nsubj/2)
                    ts = z{a};
                    ts_1 = [ts_1; ts];
                end
                for a = idx(1+(Nsubj/2):end)
                    ts = z{a};
                    ts_2 = [ts_2; ts];
                end

                [Ls1, gamma_hats1, mus1, log_likelihood1] = glmm_matlab(ts_1, iterations, k,spread,regul, delta);
                [Ls2, gamma_hats2, mus2, log_likelihood2] = glmm_matlab(ts_2, iterations, k,spread,regul, delta);

                cs = zeros(k);
                for m = 1:k
                    for n = 1:k
                        cs(m,n) = pdist2(mus1(:,m)', mus2(:,n)','cosine'); % cosine distance = 1 - cosine similarity
                    end
                end
                for m = 1:k
                    cs(m,m) = 3;
                end
                
                % Ls1, Ls2에 대해 Hungarian algorithm으로 맞는 cluster matching
                [M, ~] = matchpairs(cs, 4); %3 is cost of unassignment 

                % matching된 두 matrix의 cosine distance 중 가장 큰 매치 pair에(cosine similiary 가장 작은)
                CD_MAX(it) = max(cs(sub2ind(size(cs), M(:,1), M(:,2))));
            end
        catch
            for it = 1:iter
                idx = randperm(Nsubj);
                ts_1 = [];
                ts_2 = [];
                for a = idx(1:Nsubj/2)
                    ts = z{a};
                    ts_1 = [ts_1; ts];
                end
                for a = idx(1+(Nsubj/2):end)
                    ts = z{a};
                    ts_2 = [ts_2; ts];
                end

                [Ls1, gamma_hats1, mus1, log_likelihood1] = glmm_matlab(ts_1, iterations, k,spread,regul, delta);
                [Ls2, gamma_hats2, mus2, log_likelihood2] = glmm_matlab(ts_2, iterations, k,spread,regul, delta);

                cs = zeros(k);
                for m = 1:k
                    for n = 1:k
                        cs(m,n) = pdist2(mus1(:,m)', mus2(:,n)','cosine'); % cosine distance = 1 - cosine similarity
                    end
                end
                for m = 1:k
                    cs(m,m) = 3;
                end
                
                % Ls1, Ls2에 대해 Hungarian algorithm으로 맞는 cluster matching
                [M, ~] = matchpairs(cs, 4); %3 is cost of unassignment 

                % matching된 두 matrix의 cosine distance 중 가장 큰 매치 pair에(cosine similiary 가장 작은)
                CD_MAX(it) = max(cs(sub2ind(size(cs), M(:,1), M(:,2))));
            end
        end
        COST = mean(CD_MAX);
        cost(inD,inK) = COST;
        inD = inD+1;
        
        %save([root_path,'results/glmm/glmm_rl_k',num2str(k),'_d',num2str(delta),'.mat'], 'COST', 'Ls_iter', 'mus_iter', 'gamma_hats_iter')
    end
    inK = inK + 1;
end

%% distance(cost)가 가장 작은 hyperparameter를 선정
[d, c] = ind2sub(size(cost),find(cost == min(min(cost))));
k = K(c);
delta = Delta(d);

figure;
cost_tp = cost';
plot(cost_tp)
hold on
xticks(linspace(1, length(K), length(K)))
xticklabels(statek_name)
legend(delta_name)