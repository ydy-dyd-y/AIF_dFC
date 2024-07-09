
function L = spm_mdp_L(P,M,U,Y)
% log-likelihood function
% FORMAT L = spm_mdp_L(P,M,U,Y)
% P    - parameter structure
% M    - generative model
% U    - inputs
% Y    - observed repsonses
%
% This function runs the generative model with a given set of parameter
% values, after adding in the observations and actions on each trial
% from (real or simulated) participant data. It then sums the
% (log-)probabilities (log-likelihood) of the participant's actions under the model when it
% includes that set of parameter values. The variational Bayes fitting
% routine above uses this function to find the set of parameter values that maximize
% the probability of the participant's actions under the model (while also
% penalizing models with parameter values that move farther away from prior
% values).
%__________________________________________________________________________

if ~isstruct(P); P = spm_unvec(P,M.pE); end

% Here we re-transform parameter values out of log- or logit-space when 
% inserting them into the model to compute the log-likelihood
%--------------------------------------------------------------------------
mdp   = M.mdp;
field = fieldnames(M.pE);
for i = 1:length(field)
    if strcmp(field{i},'alpha')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'beta')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'la')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'rs')
        mdp.(field{i}) = exp(P.(field{i}));
    elseif strcmp(field{i},'eta')
        mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
    elseif strcmp(field{i},'omega')
        mdp.(field{i}) = 1/(1+exp(-P.(field{i})));
    else
        mdp.(field{i}) = exp(P.(field{i}));
    end
end

% place MDP in trial structure
%{
% --------------------------------------------------------------------------
la = mdp.la_true;  % true level of loss aversion
rs = mdp.rs_true;  % true preference magnitude for winning (higher = more risk-seeking)
%}
if isfield(M.pE,'la')&&isfield(M.pE,'rs')
    mdp.C{1} = [0 -mdp.la    ;  % Null
                0 -mdp.la;  % Loss
                0 mdp.rs ]; % win
elseif isfield(M.pE,'la')
    mdp.C{1} = [0 -mdp.la ;  % Null
                0 -mdp.la ;  % Loss
                0 0 ];       % win
elseif isfield(M.pE,'rs')
    mdp.C{1} = [0  0  ;      % Null
                0  0  ;      % Loss
                0  mdp.rs ]; % win
else
    mdp.C{1} = [0  -3  ;  % Null
                0  -3 ;  % Loss
                0  3 ]; % win
end

j = 1:numel(U); % observations for each trial
n = numel(j);   % number of trials

[MDP(1:n)] = deal(mdp);  % Create MDP with number of specified trials
[MDP.o]    = deal(U{j}); % Add observations in each trial

% solve MDP and accumulate log-likelihood
%--------------------------------------------------------------------------
MDP   = spm_MDP_VB_X(MDP); % run model with possible parameter values

L     = 0; % start (log) probability of actions given the model at 0

for i = 1:numel(Y) % Get probability of true actions for each trial
    for j = 1:numel(Y{1}(2,:)) % Only get probability of the second (controllable) state factor
        
        L = L + log(MDP(i).P(:,Y{i}(2,j),j)+ eps); % sum the (log) probabilities of each action
                                                   % given a set of possible parameter values
    end
end 

clear('MDP')

fprintf('LL: %f \n',L)