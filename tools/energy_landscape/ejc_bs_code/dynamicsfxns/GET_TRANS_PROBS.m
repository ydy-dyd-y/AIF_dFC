function [transitionProbability,transitionProbabilityMatrices,numTransitions] = GET_TRANS_PROBS(partition,subjInd,start,numClusters)

% Calculate transition probabilities in a sequence of states
% partition: integer vector, sequential cluster assignments 
% subjInd: integer vector, subject index for partition
% start: if this value is 1, then the sum of TPs whose start state is same is 1. 
%        if 2, then the sum of TPs whose next state is same is 1. if 3, then the sum of all TP is 1.
% numClusters: number of states (indexed 1:numClusters)
% Return 2D (transitionProbability) or 3D (transitionProbabilityMatrices)

partition = reshape(partition,length(partition),1); %convert to row vector
nobs = max(subjInd);
if ~exist('numClusters','var')
	numClusters = length(unique(partition));
end
possible_transitions = (numClusters)*(numClusters);
numTransitions = zeros(nobs,numClusters,numClusters);

transitionProbabilityMatrices = zeros(size(numTransitions));

for N = 1:nobs
    subjMask = subjInd == N;
    subjMask = partition(subjMask)';
    for Kinitial = 1:numClusters
        for Kfinal = 1:numClusters
            numTransitions(N,Kinitial,Kfinal) = length(strfind(subjMask,[Kinitial Kfinal]));
        end
    end
    if start == 1
        transitionProbabilityMatrices(N,:,:) = squeeze(numTransitions(N,:,:)) ./ repmat(sum(squeeze(numTransitions(N,:,:)),2),[1 numClusters]);
    elseif start == 2
        transitionProbabilityMatrices(N,:,:) = squeeze(numTransitions(N,:,:)) ./ repmat(sum(squeeze(numTransitions(N,:,:)),1)',[1 numClusters])';
    else
        transitionProbabilityMatrices(N,:,:) = squeeze(numTransitions(N,:,:)) ./ repmat(repmat(sum(sum(squeeze(numTransitions(N,:,:)),1)),[1 numClusters])', [1 numClusters]);
    end
end

transitionProbabilityMatrices(isnan(transitionProbabilityMatrices)) = 0;

transitionProbability = reshape(permute(transitionProbabilityMatrices,[1 3 2]),nobs,possible_transitions);

numTransitions = reshape(permute(numTransitions,[1 3 2]),nobs,possible_transitions);