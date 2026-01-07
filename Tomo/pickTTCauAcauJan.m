function [velPick,tPick] = pickTTCauAcauJan(SAIFFT,fInp,tInp,stnDist,vRef,vPrct)
% this function uses find peaks and outputs the travel times that closely
% match the reference traveltimes

for fNo = 1:1:length(fInp)
    %Now find the peaks
    [~,peakLoc] = findpeaks(SAIFFT(:,fNo));
    %now take the biggest peak;
    if(~isempty(peakLoc))
        tPeak = tInp(1,peakLoc);
        vPeak = stnDist./tPeak;
        vDiff = abs(vPeak-vRef(fNo,1));
        vDiffPrct = vDiff./vRef(fNo,1);
        [minVal,minInd] = min(vDiffPrct);
        if(minVal<vPrct)
            tPick(fNo,1) = tInp(1,peakLoc(minInd));
            velPick(fNo,1) = stnDist/tPick(fNo,1);
        else
            velPick(fNo,1) = NaN;
            tPick(fNo,1) = NaN;
        end
    else
        velPick(fNo,1) = NaN;
        tPick(fNo,1) = NaN;
    end
end

end