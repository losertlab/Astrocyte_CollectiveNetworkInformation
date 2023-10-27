function [h,H] = symbolTimeEvolution(states,tau,cellProbs)
% for a given symbol sequence find the entropy of all states after tau time
% steps

% INPUTS: S-symbol sequence
%        tau - time step
%       cellProbs - coarsegrained probabilities of landing in ith cell
    
    stateVals = unique(states);
    numberStates = length(unique(states));
    h = cell(numberStates,1);
    H=0;
    for s = 1:numberStates
        stateLocs=find(strcmp(states,stateVals(s)))+tau;
        % need to make sure the time step doesn't exceed size of symSeq
        stateLocs = stateLocs((stateLocs<=length(states)));
        s_tau = categorical(states(stateLocs));
        [tempProbs,~]=histcounts(s_tau,'Normalization','probability');
        numEvolStates = length(tempProbs);
        
        % calculate partitioned entropy
        tempH=0;
        for mm=1:numEvolStates
            tP = tempProbs(mm);
            tempH=tempH+tP*log(tP);
        end
        tempH=-1*tempH;
        h{s,1}=tempH;

        % calculate partitioned entropy for the entire attractor
        H=H+cellProbs(2,s)*tempH;
    end

    
end