function [tMax,phMax] = pickTTPhase(SAInp,phInp,fInp,tInp)
% this function was written to pick the phase at the observed group travel
% time

for i = 1:1:length(fInp)
    % just pick the maximum for now
    [~,maxInd] = max(SAInp(:,i));
    tMax(i,1) = tInp(maxInd);
    phMax(i,1) = phInp(maxInd,i);
end
    
end

