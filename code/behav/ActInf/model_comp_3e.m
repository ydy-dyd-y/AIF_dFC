clear all
close all      % These commands clear the workspace and close any figures
clc

root_path = '/home/dayoung/dFC/RL/';

addpath(genpath([root_path, 'tools/spm12/']))
addpath(genpath('/home/dayoung/Tool_box/spm12'))

for task_seq = 2 % 1 is loss task, 2 is reward task
    for dec_param = 1 % 1: alpha,beta, 2: beta, 3: alpha

        rng('shuffle') % This sets the random number generator to produce a different
        % random sequence each time, which leads to variability in
        % repeated simulation results (you can alse set to 'default'
        % to produce the same random sequence each time)
        
        % 1. Set up model structure
        
        T = 2;
        D{1} = [1 0]';  % {'fig.A better','fig.B better'}
        
        % For the 'behavior' state factor, we can specify that the agent always
        % begins a trial in the 'start' state (i.e., before choosing)
        
        D{2} = [1 0 0]'; % {'no choose', 'choose-fig.A','choose-fig.B'}
        
        d{1} = [.25 .25]';  % {'fig.A better','fig.B better'}
        
        % For behavior beliefs, we can specify that the agent expects with
        % certainty that it will begin a trial in the 'start' state:
        
        d{2} = [1 0 0]'; % {'no choose', 'choose-fig.A','choose-fig.B'}
        
        Ns = [length(D{1}) length(D{2})]; % number of states in each state factor (2 and 3)
        
        for i = 1:Ns(2)
            
            A{1}(:,:,i) = [1 1;  % Null
                0 0;  % Loss
                0 0]; % Win
        end
        
        pWin = .7;
        
        % 'choose-fig.A'
        A{1}(:,:,2) = [0      0;     % Null
            1-pWin pWin;  % Lose
            pWin 1-pWin]; % Win
        % 'choose-fig.B'
        A{1}(:,:,3) = [0      0;     % Null
            pWin 1-pWin;  % Lose
            1-pWin pWin]; % Win
        
        for i = 1:Ns(2)
            
            A{2}(i,:,i) = [1 1];
            
        end
        
        a{1} = A{1};
        a{2} = A{2};
        
        a{1}(:,:,2) =  [0  0;  % Null
            0 .5;  % Lose
            .5 0]; % Win
        
        a{1}(:,:,3) =  [0  0;  % Null
            0 .5;  % Lose
            .5 0]; % Win        
        
        % Controlled transitions and transition beliefs : B{:,:,u} and b(:,:,u)
        %==========================================================================
        
        B{1}(:,:,1) = [1 0;  % 'Left Better' Context
            0 1]; % 'Right Better' Context
        
        
        % Move to the Start state from any other state
        B{2}(:,:,1) = [1 1 1;  % No Choose
            0 0 0;  % choose fig.A
            0 0 0]; % choose fig.B
        
        % Move to the Choose B from any other state
        B{2}(:,:,2) = [0 0 0;  % No Choose
            1 1 1;  % choose fig.A
            0 0 0]; % choose fig.B
        
        
        % Move to the Choose Left state from any other state
        B{2}(:,:,3) = [0 0 0;  % No Choose
            0 0 0;  % choose fig.A
            1 1 1]; % choose fig.B
        
        b{1} = B{1};
        b{2} = B{2};
        % Preferred outcomes: C and c
        %==========================================================================
        
        No = [size(A{1},1) size(A{2},1)]; % number of outcomes in each outcome modality
        
        C{1} = zeros(No(1),T); % Null/Losses/Wins (3)
        C{2} = zeros(No(2),T); % Observed Behaviors (3)
        
        E = [0 1/2 1/2]';
        
        la = 3; % loss aversivity
        rs = 3; % reward sensitivity
        
        if task_seq == 1  % loss task
            C{1}(:,:) =    [0 -la;  % Null  %unleasness occurred by not getting more information is calculated in other parts
                0 -la;  % Loss
                0 0]; % win
        else
            C{1}(:,:) =    [0 0;  % Null
                0 0  ;  % Loss
                0 rs]; % win
        end
        
        Np = 3; % Number of policies
        Nf = 2; % Number of state factors
        
        V         = ones(T-1,Np,Nf);
        
        V(:,:,1) = [1 1 1]; % Context state is not controllable
        
        V(:,:,2) = [1 2 3];     % action
        
        
        
        % Additional optional parameters.
        %==========================================================================
        
        eta = 0.5; % By default we here set this to 0.5, but try changing its value
        % to see how it affects model behavior
        
        
        omega = 1; % By default we here set this to 1 (indicating no forgetting,
        
        beta = 1; % By default this is set to 1, but try increasing its value
        
        alpha = 32;  % Any positive number. 1 is very low, 32 is fairly high;
        
        
        erp = 1; % By default we here set this to 1, but try increasing its value
        
        tau = 12; % Here we set this to 12 to simulate smooth physiological responses,
        % but try adjusting its value to see how it affects simulated
        % neural (and behavioral) responses
        
        
        % 2. Define MDP Structure
        %==========================================================================
        %==========================================================================
        % 1. including d
        mdp1.T = T;                    % Number of time steps
        mdp1.V = V;                    % allowable (deep) policies
        
        %mdp.U = U;                   % We could have instead used shallow
        % policies (specifying U instead of V).
        
        mdp1.A = A;                    % state-outcome mapping
        mdp1.B = B;                    % transition probabilities
        mdp1.C = C;                    % preferred states
        mdp1.D = D;                    % priors over initial states
        mdp1.E = E;
        
        mdp1.d = d;                    % enable learning priors over initial states
        mdp1.un = [1, 0; 0, 0.5; 0, 0.5];
        
        mdp1.eta = eta;                % learning rate
        mdp1.omega = omega;            % forgetting rate
        mdp1.alpha = alpha;            % action precision
        mdp1.beta = beta;              % expected precision of expected free energy over policies
        mdp1.erp = erp;                % degree of belief resetting at each timestep
        mdp1.tau = tau;                % time constant for evidence accumulation
        
        label.factor{1}   = 'contexts';   label.name{1}    = {'fig.A-better','fig.B-better'};
        label.factor{2}   = 'choice states';     label.name{2}    = {'No choose','choose fig.A','choose fig.B'};
        label.modality{1} = 'win/lose';  label.outcome{1} = {'null','lose','win'};
        label.modality{2} = 'observed action';  label.outcome{2} = {'No choose','choose fig.A','choose fig.B'};
        label.action{2} = {'No choose','choose fig.A','choose fig.B'};
        mdp1.label = label;
        mdp1 = spm_MDP_check(mdp1);
        
        % 3. including a, d
        mdp3.T = T;                    % Number of time steps
        mdp3.V = V;                    % allowable (deep) policies
        
        %mdp.U = U;                   % We could have instead used shallow
        % policies (specifying U instead of V).
        
        mdp3.A = A;                    % state-outcome mapping
        mdp3.B = B;                    % transition probabilities
        mdp3.C = C;                    % preferred states
        mdp3.D = D;                    % priors over initial states
        mdp3.E = E;
        
        mdp3.a = a;
        mdp3.d = d;                    % enable learning priors over initial states
        mdp3.un = [1, 0; 0, 0.5; 0, 0.5];
        
        mdp3.eta = eta;                % learning rate
        mdp3.omega = omega;            % forgetting rate
        mdp3.alpha = alpha;            % action precision
        mdp3.beta = beta;              % expected precision of expected free energy over policies
        mdp3.erp = erp;                % degree of belief resetting at each timestep
        mdp3.tau = tau;                % time constant for evidence accumulation
        
        label.factor{1}   = 'contexts';   label.name{1}    = {'fig.A-better','fig.B-better'};
        label.factor{2}   = 'choice states';     label.name{2}    = {'No choose','choose fig.A','choose fig.B'};
        label.modality{1} = 'win/lose';  label.outcome{1} = {'null','lose','win'};
        label.modality{2} = 'observed action';  label.outcome{2} = {'No choose','choose fig.A','choose fig.B'};
        label.action{2} = {'No choose','choose fig.A','choose fig.B'};
        mdp3.label = label;
        mdp3 = spm_MDP_check(mdp3);
        
        clear beta
        clear alpha
        clear eta
        clear omega
        clear la
        clear rs
        
        load([root_path, 'data/behav/behavior_data.mat'])
        
        %% 7. Model comparison
        %==========================================================================
        
        % Create vectors/matrices that will store results
        
        F_3e_L1_params = [];
        F_3e_L3_params = [];
        
        avg_LL_3e_L1_params = [];
        avg_prob_3e_L1_params = [];
        avg_LL_3e_L3_params = [];
        avg_prob_3e_L3_params = [];
        
        GCM_3e_L1 = {};
        GCM_3e_L3 = {};
        
        GMDP_3e_L1 = {};
        GMDP_3e_L3 = {};
        
        Sim_params_3e_L1 = [];
        Sim_params_3e_L3 = [];
        
        % Set up reversal learning trials like before
        
        N = 40; % number of trials
        MDP = mdp1;
        [MDP(1:N)] = deal(MDP);
        
        ncores=8;
        %{
        try
            pool = parpool('local',ncores);
            parfor si = 1:length(data)  % specify different true risk-seeking values (prior = 2)
                MDP_temp = MDP;
                for ti = 1:N
                    MDP_temp(ti).u = ones(Nf,1);  % u : [Num_state factor, performing action]    % action -> feed -> fix
                    MDP_temp(ti).o = ones(Nf,T);    % o : [Num_outcome modality, T]  % cue -> action -> feed -> fix
                    if data{1,si}.response{1,task_seq}(ti) == 0
                        MDP_temp(ti).o(2,2) = 3;
                        MDP_temp(ti).u(2,1) = 3;    % choose non-better choice(Fig.B)
                    else
                        MDP_temp(ti).o(2,2) = 2;
                        MDP_temp(ti).u(2,1) = 2;    % choose better choice(Fig.A)
                    end
                    
                    if data{1,si}.input{1,task_seq}(ti) == data{1,si}.response{1,task_seq}(ti)
                        MDP_temp(ti).o(1,2) = 3;   % win
                    else
                        MDP_temp(ti).o(1,2) = 2;   % lose
                    end
                end
                
                %MDP_temp = spm_MDP_VB_X(MDP_temp);
                
                DCM1.MDP   = mdp1;                  % MDP model that will be estimated
                if task_seq == 1
                    if dec_param == 1
                        DCM1.field = {'alpha', 'beta', 'la', 'eta'}; % parameter (field) names to optimise
                    elseif dec_param == 2
                        DCM1.field = {'beta', 'la', 'eta'};
                    else
                        DCM1.field = {'alpha', 'la', 'eta'};
                    end
                else
                    if dec_param == 1
                        DCM1.field = {'alpha', 'beta', 'rs', 'eta'}; % parameter (field) names to optimise
                    elseif dec_param == 2
                        DCM1.field = {'beta', 'rs', 'eta'};
                    else
                        DCM1.field = {'alpha', 'rs', 'eta'};
                    end
                end
                
                DCM1.U     = {MDP_temp.o};              % include the observations made by (real
                % or simulated) participants
                
                DCM1.Y     = {MDP_temp.u};              % include the actions made by (real or
                % simulated) participants
                
                DCM1       = Estimate_parameters(DCM1); % Run the parameter estimation function
                
                % Convert parameters back out of log- or logit-space
                
                field = fieldnames(DCM1.M.pE);
                for i = 1:length(field)
                    if strcmp(field{i},'eta')
                        DCM1.prior(i) = 1/(1+exp(-DCM1.M.pE.(field{i})));
                        DCM1.posterior(i) = 1/(1+exp(-DCM1.Ep.(field{i})));
                    elseif strcmp(field{i},'omega')
                        DCM1.prior(i) = 1/(1+exp(-DCM1.M.pE.(field{i})));
                        DCM1.posterior(i) = 1/(1+exp(-DCM1.Ep.(field{i})));
                    else
                        DCM1.prior(i) = exp(DCM1.M.pE.(field{i}));
                        DCM1.posterior(i) = exp(DCM1.Ep.(field{i}));
                    end
                end
                
                F_3e_L1_params = [F_3e_L1_params DCM1.F];% Get free energies for each participant's model
                
                GCM_3e_L1   = [GCM_3e_L1;{DCM1}]; % Save DCM for each participant
                
                % Get Log-likelihood and action probabilities for best-fit model
                
                MDP_best = MDP;
                
                if dec_param == 1
                    [MDP_best(1:N).alpha] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).beta] = deal(DCM1.posterior(2));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(4));
                elseif dec_param == 2
                    [MDP_best(1:N).beta] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(3));
                else
                    [MDP_best(1:N).alpha] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(3));
                end
                
                if task_seq == 1
                    if dec_param == 1
                        C_fit_best = [0 -DCM1.posterior(3)   ;  % Null
                            0 -DCM1.posterior(3);  % Loss
                            0 0 ]; % win
                    else
                        C_fit_best = [0 -DCM1.posterior(2)   ;  % Null
                            0 -DCM1.posterior(2);  % Loss
                            0 0 ]; % win
                    end
                else
                    if dec_param == 1
                        C_fit_best = [0 0   ;  % Null
                            0 0   ;  % Loss
                            0 DCM1.posterior(3) ]; % win
                    else
                        C_fit_best = [0 0  ;  % Null
                            0 0  ;  % Loss
                            0 DCM1.posterior(2) ]; % win
                    end
                end
                
                for i = 1:N
                    MDP_best(i).C{1} = C_fit_best;
                end
                
                for i = 1:N
                    MDP_best(i).o = MDP_temp(i).o;
                end
                
                for i = 1:N
                    MDP_best(i).u = MDP_temp(i).u;
                end
                
                MDP_best   = spm_MDP_VB_X(MDP_best); % run model with best parameter values
                
                % Get sum of log-likelihoods for each action across trials
                
                L     = 0; % start (log) probability of actions given the model at 0
                total_prob = 0;
                
                for i = 1:numel(MDP_best) % Get probability of true actions for each trial
                    for j = 1:numel(MDP_best(1).u(2,:)) % Only get probability of the second (controllable) state factor
                        
                        L = L + log(MDP_best(i).P(:,MDP_best(i).u(2,j),j)+ eps); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        total_prob = total_prob + MDP_best(i).P(:,MDP_best(i).u(2,j),j); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        
                    end
                end
                
                % Get the average log-likelihood for each participant and average action
                % probability of each participant under best-fit parameters
                
                avg_LL_3e_L1 = L/(size(MDP_best,2)*2);
                
                avg_LL_3e_L1_params = [avg_LL_3e_L1_params; avg_LL_3e_L1];
                
                avg_prob_3e_L1 = total_prob/(size(MDP_best,2)*2);
                
                avg_prob_3e_L1_params = [avg_prob_3e_L1_params; avg_prob_3e_L1];
                
                % Store true and estimated parameters to assess recoverability
                
                Sim_params_3e_L1 = [Sim_params_3e_L1; DCM1.posterior];% Get posteriors
                
                GMDP_3e_L1 = [GMDP_3e_L1;{MDP_best}];
                
            end
        catch
            disp('Parallel Computing Toolbox is not available')
            for si = 1:length(data)  % specify different true risk-seeking values (prior = 2)
                
                MDP_temp = MDP;
                for ti = 1:N
                    MDP_temp(ti).u = ones(Nf,1);  % u : [Num_state factor, performing action]    % action -> feed -> fix
                    MDP_temp(ti).o = ones(Nf,T);    % o : [Num_outcome modality, T]  % cue -> action -> feed -> fix
                    if data{1,si}.response{1,task_seq}(ti) == 0
                        MDP_temp(ti).o(2,2) = 3;
                        MDP_temp(ti).u(2,1) = 3;    % choose non-better choice(Fig.B)
                    else
                        MDP_temp(ti).o(2,2) = 2;
                        MDP_temp(ti).u(2,1) = 2;    % choose better choice(Fig.A)
                    end
                    
                    if data{1,si}.input{1,task_seq}(ti) == data{1,si}.response{1,task_seq}(ti)
                        MDP_temp(ti).o(1,2) = 3;   % win
                    else
                        MDP_temp(ti).o(1,2) = 2;   % lose
                    end
                end
                
                %MDP_temp = spm_MDP_VB_X(MDP_temp);
                
                DCM1.MDP   = mdp1;                  % MDP model that will be estimated
                
                if task_seq == 1
                    if dec_param == 1
                        DCM1.field = {'alpha', 'beta', 'la', 'eta'}; % parameter (field) names to optimise
                    elseif dec_param == 2
                        DCM1.field = {'beta', 'la', 'eta'};
                    else
                        DCM1.field = {'alpha', 'la', 'eta'};
                    end
                else
                    if dec_param == 1
                        DCM1.field = {'alpha', 'beta', 'rs', 'eta'}; % parameter (field) names to optimise
                    elseif dec_param == 2
                        DCM1.field = {'beta', 'rs', 'eta'};
                    else
                        DCM1.field = {'alpha', 'rs', 'eta'};
                    end
                end
                
                DCM1.U     = {MDP_temp.o};              % include the observations made by (real
                % or simulated) participants
                
                DCM1.Y     = {MDP_temp.u};              % include the actions made by (real or
                % simulated) participants
                
                DCM1       = Estimate_parameters(DCM1); % Run the parameter estimation function
                
                % Convert parameters back out of log- or logit-space
                
                field = fieldnames(DCM1.M.pE);
                for i = 1:length(field)
                    if strcmp(field{i},'eta')
                        DCM1.prior(i) = 1/(1+exp(-DCM1.M.pE.(field{i})));
                        DCM1.posterior(i) = 1/(1+exp(-DCM1.Ep.(field{i})));
                    elseif strcmp(field{i},'omega')
                        DCM1.prior(i) = 1/(1+exp(-DCM1.M.pE.(field{i})));
                        DCM1.posterior(i) = 1/(1+exp(-DCM1.Ep.(field{i})));
                    else
                        DCM1.prior(i) = exp(DCM1.M.pE.(field{i}));
                        DCM1.posterior(i) = exp(DCM1.Ep.(field{i}));
                    end
                end
                
                F_3e_L1_params = [F_3e_L1_params DCM1.F];% Get free energies for each participant's model
                
                GCM_3e_L1   = [GCM_3e_L1;{DCM1}]; % Save DCM for each participant
                
                % Get Log-likelihood and action probabilities for best-fit model
                
                MDP_best = MDP;
                
                if dec_param == 1
                    [MDP_best(1:N).alpha] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).beta] = deal(DCM1.posterior(2));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(4));
                elseif dec_param == 2
                    [MDP_best(1:N).beta] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(3));
                else
                    [MDP_best(1:N).alpha] = deal(DCM1.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM1.posterior(3));
                end
                
                if task_seq == 1
                    if dec_param == 1
                        C_fit_best = [0 -DCM1.posterior(3)   ;  % Null
                            0 -DCM1.posterior(3);  % Loss
                            0 0 ]; % win
                    else
                        C_fit_best = [0 -DCM1.posterior(2)   ;  % Null
                            0 -DCM1.posterior(2);  % Loss
                            0 0 ]; % win
                    end
                else
                    if dec_param == 1
                        C_fit_best = [0 0   ;  % Null
                            0 0   ;  % Loss
                            0 DCM1.posterior(3) ]; % win
                    else
                        C_fit_best = [0 0  ;  % Null
                            0 0  ;  % Loss
                            0 DCM1.posterior(2) ]; % win
                    end
                end
                
                for i = 1:N
                    MDP_best(i).C{1} = C_fit_best;
                end
                
                for i = 1:N
                    MDP_best(i).o = MDP_temp(i).o;
                end
                
                for i = 1:N
                    MDP_best(i).u = MDP_temp(i).u;
                end
                
                MDP_best   = spm_MDP_VB_X(MDP_best); % run model with best parameter values
                
                % Get sum of log-likelihoods for each action across trials
                
                L     = 0; % start (log) probability of actions given the model at 0
                total_prob = 0;
                
                for i = 1:numel(MDP_best) % Get probability of true actions for each trial
                    for j = 1:numel(MDP_best(1).u(2,:)) % Only get probability of the second (controllable) state factor
                        
                        L = L + log(MDP_best(i).P(:,MDP_best(i).u(2,j),j)+ eps); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        total_prob = total_prob + MDP_best(i).P(:,MDP_best(i).u(2,j),j); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        
                    end
                end
                
                % Get the average log-likelihood for each participant and average action
                % probability of each participant under best-fit parameters
                
                avg_LL_3e_L1 = L/(size(MDP_best,2)*2);
                
                avg_LL_3e_L1_params = [avg_LL_3e_L1_params; avg_LL_3e_L1];
                
                avg_prob_3e_L1 = total_prob/(size(MDP_best,2)*2);
                
                avg_prob_3e_L1_params = [avg_prob_3e_L1_params; avg_prob_3e_L1];
                
                % Store true and estimated parameters to assess recoverability
                
                Sim_params_3e_L1 = [Sim_params_3e_L1; DCM1.posterior];% Get posteriors
                
                GMDP_3e_L1 = [GMDP_3e_L1;{MDP_best}];
                
            end
        end
        
        % Separately store true and simulated parameters
        if dec_param == 1
            Estimated_alpha_3e_L1 = Sim_params_3e_L1(:,1);
            Estimated_beta_3e_L1 = Sim_params_3e_L1(:,2);
            if task_seq == 1
                Estimated_la_3e_L1 = Sim_params_3e_L1(:,3);
            else
                Estimated_rs_3e_L1 = Sim_params_3e_L1(:,3);
            end
            Estimated_eta_3e_L1 = Sim_params_3e_L1(:,4);
        elseif dec_param == 2
            Estimated_beta_3e_L1 = Sim_params_3e_L1(:,1);
            if task_seq == 1
                Estimated_la_3e_L1 = Sim_params_3e_L1(:,2);
            else
                Estimated_rs_3e_L1 = Sim_params_3e_L1(:,2);
            end
            Estimated_eta_3e_L1 = Sim_params_3e_L1(:,3);
        else
            Estimated_alpha_3e_L1 = Sim_params_3e_L1(:,1);
            if task_seq == 1
                Estimated_la_3e_L1 = Sim_params_3e_L1(:,2);
            else
                Estimated_rs_3e_L1 = Sim_params_3e_L1(:,2);
            end
            Estimated_eta_3e_L1 = Sim_params_3e_L1(:,3);
        end
        %}

        %% Generate free energies for model fits for L3
        MDP = mdp3;
        [MDP(1:N)] = deal(MDP);
        
        try
            pool = parpool('local',ncores);
        end         
        tmp = exist('pool');
        if tmp == 1
            parfor si = 1:length(data)  % specify different true risk-seeking values (prior = 2)
                
                MDP_temp = MDP;
                
                for ti = 1:N
                    MDP_temp(ti).u = ones(Nf,1);  % u : [Num_state factor, performing action]    % action -> feed -> fix
                    MDP_temp(ti).o = ones(Nf,T);    % o : [Num_outcome modality, T]  % cue -> action -> feed -> fix
                    if data{1,si}.response{1,task_seq}(ti) == 0
                        MDP_temp(ti).o(2,2) = 3;
                        MDP_temp(ti).u(2,1) = 3;    % choose non-better choice(Fig.B)
                    else
                        MDP_temp(ti).o(2,2) = 2;
                        MDP_temp(ti).u(2,1) = 2;    % choose better choice(Fig.A)
                    end
                    
                    if data{1,si}.input{1,task_seq}(ti) == data{1,si}.response{1,task_seq}(ti)
                        MDP_temp(ti).o(1,2) = 3;   % win
                    else
                        MDP_temp(ti).o(1,2) = 2;   % lose
                    end
                end
                
                %MDP_temp = spm_MDP_VB_X(MDP_temp);
                
                DCM3.MDP   = mdp3;                  % MDP model that will be estimated
                
                if task_seq == 1
                    if dec_param == 1
                        DCM3.field = {'alpha', 'beta', 'la', 'eta'};
                    elseif dec_param == 2
                        DCM3.field = {'beta', 'la', 'eta'};
                    else
                        DCM3.field = {'alpha', 'la', 'eta'};
                    end
                else
                    if dec_param == 1
                        DCM3.field = {'alpha', 'beta', 'rs', 'eta'};
                    elseif dec_param == 2
                        DCM3.field = {'beta', 'rs', 'eta'};
                    else
                        DCM3.field = {'alpha', 'rs', 'eta'};
                    end
                end
                
                DCM3.U     = {MDP_temp.o};        % include the observations made by (real
                % or simulated) participants
                
                DCM3.Y     = {MDP_temp.u};        % include the actions made by (real or
                % simulated) participants
                
                DCM3       = Estimate_parameters(DCM3); % Run the parameter estimation function
                
                % Convert parameters back out of log- or logit-space
                
                field = fieldnames(DCM3.M.pE);
                for i = 1:length(field)
                    if strcmp(field{i},'eta')
                        DCM3.prior(i) = 1/(1+exp(-DCM3.M.pE.(field{i})));
                        DCM3.posterior(i) = 1/(1+exp(-DCM3.Ep.(field{i})));
                    elseif strcmp(field{i},'omega')
                        DCM3.prior(i) = 1/(1+exp(-DCM3.M.pE.(field{i})));
                        DCM3.posterior(i) = 1/(1+exp(-DCM3.Ep.(field{i})));
                    else
                        DCM3.prior(i) = exp(DCM3.M.pE.(field{i}));
                        DCM3.posterior(i) = exp(DCM3.Ep.(field{i}));
                    end
                end
                
                F_3e_L3_params = [F_3e_L3_params DCM3.F]; % Get free energies for each participant's model
                
                GCM_3e_L3   = [GCM_3e_L3;{DCM3}]; % Save DCM for each participant
                
                % Get Log-likelihood and action probabilities for best-fit model
                
                MDP_best = MDP;
                if dec_param == 1
                    [MDP_best(1:N).alpha] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).beta] = deal(DCM3.posterior(2));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(4));
                elseif dec_param == 2
                    [MDP_best(1:N).beta] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(3));
                else
                    [MDP_best(1:N).alpha] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(3));
                end
                
                if task_seq == 1
                    if dec_param == 1
                        C_fit_best = [0 -DCM3.posterior(3)   ;  % Null
                            0 -DCM3.posterior(3);  % Loss
                            0 0 ]; % win
                    else
                        C_fit_best = [0 -DCM3.posterior(2)   ;  % Null
                            0 -DCM3.posterior(2);  % Loss
                            0 0 ]; % win
                    end
                else
                    if dec_param == 1
                        C_fit_best = [0 0   ;  % Null
                            0 0   ;  % Loss
                            0 DCM3.posterior(3) ]; % win
                    else
                        C_fit_best = [0 0  ;  % Null
                            0 0  ;  % Loss
                            0 DCM3.posterior(2) ]; % win
                    end
                end
                
                for i = 1:N
                    MDP_best(i).C{1} = C_fit_best;
                end
                
                for i = 1:N
                    MDP_best(i).o = MDP_temp(i).o;
                end
                
                for i = 1:N
                    MDP_best(i).u = MDP_temp(i).u;
                end
                
                MDP_best   = spm_MDP_VB_X(MDP_best); % run model with best parameter values
                
                % Get sum of log-likelihoods for each action across trials
                
                L     = 0; % start (log) probability of actions given the model at 0
                total_prob = 0;
                
                for i = 1:numel(MDP_best) % Get probability of true actions for each trial
                    for j = 1:numel(MDP_best(1).u(2,:)) % Only get probability of the second (controllable) state factor
                        
                        L = L + log(MDP_best(i).P(:,MDP_best(i).u(2,j),j)+ eps); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        total_prob = total_prob + MDP_best(i).P(:,MDP_best(i).u(2,j),j); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        
                    end
                end
                
                % Get the average log-likelihood for each participant and average action
                % probability of each participant under best-fit parameters
                
                avg_LL_3e_L3 = L/(size(MDP_best,2)*2);
                avg_LL_3e_L3_params = [avg_LL_3e_L3_params; avg_LL_3e_L3];
                avg_prob_3e_L3 = total_prob/(size(MDP_best,2)*2);
                avg_prob_3e_L3_params = [avg_prob_3e_L3_params; avg_prob_3e_L3];
                
                % Store true and estimated parameters to assess recoverability
                
                Sim_params_3e_L3 = [Sim_params_3e_L3; DCM3.posterior];% Get posteriors
                GMDP_3e_L3 = [GMDP_3e_L3;{MDP_best}];
                
            end
        else
            for si = 1:length(data)  % specify different true risk-seeking values (prior = 2)
                
                MDP_temp = MDP;
                
                for ti = 1:N
                    MDP_temp(ti).u = ones(Nf,1);  % u : [Num_state factor, performing action]    % action -> feed -> fix
                    MDP_temp(ti).o = ones(Nf,T);    % o : [Num_outcome modality, T]  % cue -> action -> feed -> fix
                    if data{1,si}.response{1,task_seq}(ti) == 0
                        MDP_temp(ti).o(2,2) = 3;
                        MDP_temp(ti).u(2,1) = 3;    % choose non-better choice(Fig.B)
                    else
                        MDP_temp(ti).o(2,2) = 2;
                        MDP_temp(ti).u(2,1) = 2;    % choose better choice(Fig.A)
                    end
                    
                    if data{1,si}.input{1,task_seq}(ti) == data{1,si}.response{1,task_seq}(ti)
                        MDP_temp(ti).o(1,2) = 3;   % win
                    else
                        MDP_temp(ti).o(1,2) = 2;   % lose
                    end
                end
                
                %MDP_temp = spm_MDP_VB_X(MDP_temp);
                
                DCM3.MDP   = mdp3;                  % MDP model that will be estimated
                
                if task_seq == 1
                    if dec_param == 1
                        DCM3.field = {'alpha', 'beta', 'la', 'eta'};
                    elseif dec_param == 2
                        DCM3.field = {'beta', 'la', 'eta'};
                    else
                        DCM3.field = {'alpha', 'la', 'eta'};
                    end
                else
                    if dec_param == 1
                        DCM3.field = {'alpha', 'beta', 'rs', 'eta'};
                    elseif dec_param == 2
                        DCM3.field = {'beta', 'rs', 'eta'};
                    else
                        DCM3.field = {'alpha', 'rs', 'eta'};
                    end
                end
                
                DCM3.U     = {MDP_temp.o};        % include the observations made by (real
                % or simulated) participants
                
                DCM3.Y     = {MDP_temp.u};        % include the actions made by (real or
                % simulated) participants
                
                DCM3       = Estimate_parameters(DCM3); % Run the parameter estimation function
                
                % Convert parameters back out of log- or logit-space
                
                field = fieldnames(DCM3.M.pE);
                for i = 1:length(field)
                    if strcmp(field{i},'eta')
                        DCM3.prior(i) = 1/(1+exp(-DCM3.M.pE.(field{i})));
                        DCM3.posterior(i) = 1/(1+exp(-DCM3.Ep.(field{i})));
                    elseif strcmp(field{i},'omega')
                        DCM3.prior(i) = 1/(1+exp(-DCM3.M.pE.(field{i})));
                        DCM3.posterior(i) = 1/(1+exp(-DCM3.Ep.(field{i})));
                    else
                        DCM3.prior(i) = exp(DCM3.M.pE.(field{i}));
                        DCM3.posterior(i) = exp(DCM3.Ep.(field{i}));
                    end
                end
                
                F_3e_L3_params = [F_3e_L3_params DCM3.F]; % Get free energies for each participant's model
                
                GCM_3e_L3   = [GCM_3e_L3;{DCM3}]; % Save DCM for each participant
                
                % Get Log-likelihood and action probabilities for best-fit model
                
                MDP_best = MDP;
                if dec_param == 1
                    [MDP_best(1:N).alpha] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).beta] = deal(DCM3.posterior(2));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(4));
                elseif dec_param == 2
                    [MDP_best(1:N).beta] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(3));
                else
                    [MDP_best(1:N).alpha] = deal(DCM3.posterior(1));
                    [MDP_best(1:N).eta] = deal(DCM3.posterior(3));
                end
                
                if task_seq == 1
                    if dec_param == 1
                        C_fit_best = [0 -DCM3.posterior(3)   ;  % Null
                            0 -DCM3.posterior(3);  % Loss
                            0 0 ]; % win
                    else
                        C_fit_best = [0 -DCM3.posterior(2)   ;  % Null
                            0 -DCM3.posterior(2);  % Loss
                            0 0 ]; % win
                    end
                else
                    if dec_param == 1
                        C_fit_best = [0 0   ;  % Null
                            0 0   ;  % Loss
                            0 DCM3.posterior(3) ]; % win
                    else
                        C_fit_best = [0 0  ;  % Null
                            0 0  ;  % Loss
                            0 DCM3.posterior(2) ]; % win
                    end
                end
                
                for i = 1:N
                    MDP_best(i).C{1} = C_fit_best;
                end
                
                for i = 1:N
                    MDP_best(i).o = MDP_temp(i).o;
                end
                
                for i = 1:N
                    MDP_best(i).u = MDP_temp(i).u;
                end
                
                MDP_best   = spm_MDP_VB_X(MDP_best); % run model with best parameter values
                
                % Get sum of log-likelihoods for each action across trials
                
                L     = 0; % start (log) probability of actions given the model at 0
                total_prob = 0;
                
                for i = 1:numel(MDP_best) % Get probability of true actions for each trial
                    for j = 1:numel(MDP_best(1).u(2,:)) % Only get probability of the second (controllable) state factor
                        
                        L = L + log(MDP_best(i).P(:,MDP_best(i).u(2,j),j)+ eps); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        total_prob = total_prob + MDP_best(i).P(:,MDP_best(i).u(2,j),j); % sum the (log) probabilities of each action
                        % given a set of possible parameter values
                        
                    end
                end
                
                % Get the average log-likelihood for each participant and average action
                % probability of each participant under best-fit parameters
                
                avg_LL_3e_L3 = L/(size(MDP_best,2)*2);
                avg_LL_3e_L3_params = [avg_LL_3e_L3_params; avg_LL_3e_L3];
                avg_prob_3e_L3 = total_prob/(size(MDP_best,2)*2);
                avg_prob_3e_L3_params = [avg_prob_3e_L3_params; avg_prob_3e_L3];
                
                % Store true and estimated parameters to assess recoverability
                
                Sim_params_3e_L3 = [Sim_params_3e_L3; DCM3.posterior];% Get posteriors
                GMDP_3e_L3 = [GMDP_3e_L3;{MDP_best}];
                
            end
        end
        % Separately store true and simulated parameters
        if dec_param == 1
            Estimated_alpha_3e_L3 = Sim_params_3e_L3(:,1);
            Estimated_beta_3e_L3 = Sim_params_3e_L3(:,2);
            if task_seq == 1
                Estimated_la_3e_L3 = Sim_params_3e_L3(:,3);
            else
                Estimated_rs_3e_L3 = Sim_params_3e_L3(:,3);
            end
            Estimated_eta_3e_L3 = Sim_params_3e_L3(:,4);
        elseif dec_param == 2
            Estimated_beta_3e_L3 = Sim_params_3e_L3(:,1);
            if task_seq == 1
                Estimated_la_3e_L3 = Sim_params_3e_L3(:,2);
            else
                Estimated_rs_3e_L3 = Sim_params_3e_L3(:,2);
            end
            Estimated_eta_3e_L3 = Sim_params_3e_L3(:,3);
        else
            Estimated_alpha_3e_L3 = Sim_params_3e_L3(:,1);
            if task_seq == 1
                Estimated_la_3e_L3 = Sim_params_3e_L3(:,2);
            else
                Estimated_rs_3e_L3 = Sim_params_3e_L3(:,2);
            end
            Estimated_eta_3e_L3 = Sim_params_3e_L3(:,3);
        end
        
        %%
        clear alpha
        
        % Random Effects Bayesian Model Comparison (of Free Energies of best-fit
        % models per participant):
        
        F_3e_L3_params = F_3e_L3_params'; % Convert free energies to column vectors
        
        average_LL_3e_L3 = mean(avg_LL_3e_L3_params);
        average_action_probability_3e_L3 = mean(avg_prob_3e_L3_params);
                
        if dec_param == 1
            if task_seq == 1
                save(['/u6/pmdd/behav/loss_est_3e.mat'], 'Estimated_alpha_3e_L3', ...
                    'Estimated_beta_3e_L3', ...
                    'Estimated_la_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/loss_result_3e.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            else
                save(['/u6/pmdd/behav/reward_est_3e.mat'], 'Estimated_alpha_3e_L3', ...
                    'Estimated_beta_3e_L3', ...
                    'Estimated_rs_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/reward_result_3e.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            end
        elseif dec_param == 2
            if task_seq == 1
                save(['/u6/pmdd/behav/loss_est_3e_be.mat'], 'Estimated_beta_3e_L3', ...
                    'Estimated_la_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/loss_result_3e_be.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            else
                save(['/u6/pmdd/behav/reward_est_3e_be.mat'], 'Estimated_beta_3e_L3', ...
                    'Estimated_rs_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/reward_result_3e_be.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            end
        else
            if task_seq == 1
                save(['/u6/pmdd/behav/loss_est_3e_al.mat'], 'Estimated_alpha_3e_L3', ...
                    'Estimated_la_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/loss_result_3e_al.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            else
                save(['/u6/pmdd/behav/reward_est_3e_al.mat'], 'Estimated_alpha_3e_L3', ...
                    'Estimated_rs_3e_L3', 'Estimated_eta_3e_L3')
                
                save(['/u6/pmdd/behav/reward_result_3e_al.mat'], 'F_3e_L3_params', 'GCM_3e_L3', ...
                    'GMDP_3e_L3', 'average_action_probability_3e_L3', ...
                    'average_LL_3e_L3', 'avg_LL_3e_L3_params', 'avg_prob_3e_L3_params')
            end
        end
        clearvars('-except', 'task_seq', 'root_path');
        %close all      % These commands clear the workspace and close any figures
        %clc
    end
end