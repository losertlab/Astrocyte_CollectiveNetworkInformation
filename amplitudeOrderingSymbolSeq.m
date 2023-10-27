function [S,cellCounts,probStates] = amplitudeOrderingSymbolSeq(trace,m,l)
% using ampltiude sorting for a m dim with time delay l construct a symbol sequence
% from x(t) time series. m! number of permutations
    numFrames = length(trace);
    S = [];


    % last potential index
    rem = mod((numFrames-1),l);
    lpi = numFrames-rem;
    while lpi+(m-1)>numFrames
        lpi=lpi-(m-1);
    end
    timeSteps = 1:l:lpi;
    numSymb = length(timeSteps);
    S=cell(numSymb,1);
    idx=1;
    for ii=timeSteps
        segTrace = trace(ii:ii+(m-1));
        [~,I]=sort(segTrace);
        sym=[];
        for kk=1:m
            sym = strcat(sym,int2str(I(kk)));
        end
        S{idx}=sym;
        idx=idx+1;
    end


    % idx=1; %brute force method
    % while l*(idx-1)+1<numFrames-m
    %     tt = l*(idx-1)+1;
    %     segTrace = trace(tt:tt+(m-1));
    %     [~,I]=sort(segTrace);
    %     sym=[];
    %     for kk=1:m
    %         sym = strcat(sym,int2str(I(kk)));
    %     end
    %     S{idx}=sym;
    %     idx=idx+1;
    % end
    S=convertCharsToStrings(S)';
    % get probabilities
    symSeq = categorical(S);
    % histogram(ttwoStateCat,'Normalization','probability');
    [cellCounts,probStates]=histcounts(symSeq);
    cellCounts(2,:)=cellCounts/numSymb;
    % cellProbs is the coarse grained probability of the system falling
    % into the ith cell... Consult:
    % Kota Shiozawa, Taisuke Uemura, and Isao T. Tokuda, Phys. Rev. E 107, 01420
    % Cell here does not refer to biological cell but a phase space

end