function [tPick,tErrOut] = pickTT(SAInp,fIn,tIn,thTime,thFreq)
% this function was written to extract the travel times after FTAN
% there is ofcourse a guiding travel-time to ensure reliable picks

% interpolate the theoretical times on fIn
useT = interp1(thFreq,thTime,fIn);

for i = 1:1:length(fIn)
    [pks,locs] = findpeaks(SAInp(:,i));
    tError = [];
    tNow = [];
    for k = 1:1:length(locs)
        tNow(k,1) = tIn(locs(k,1));
        tDiff = abs(tNow(k,1)-useT(i,1));
        tError(k,1) = tDiff/useT(i,1)*100;
    end
    
    [goodInd] = find(tError<=20);
    if(~isempty(goodInd))
        tNowGood = tNow(goodInd,1);
        [~,minInd] = min(tError(goodInd,1));
        tPick(i,1) = tNowGood(minInd,1);
        tErrOut(i,1) = tError(minInd,1);
    else
        tPick(i,1) = NaN;
        tErrOut(i,1) = NaN;
    end 
end

