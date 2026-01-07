function [velPickUpdated,timePickUpdated] = pickTTSymmJan(SAIFFT,fInp,tInp,stnDist,vBand,fSamp)
% this function applies findpeaks to find the peaks in the FTAN image
% give out the velocities for the symmetric cross-correlation

% generate a tukey window to taper the end of SAIFFT
tkySmall = tukeywin(100,1);
winVal = [ones(451,1);tkySmall(51:100,1)];

tLow = stnDist./vBand(:,2);
tHigh = stnDist./vBand(:,3);

tLow = round(tLow*fSamp);
tHigh = round(tHigh*fSamp);

tLow(tLow==0) = 1;
tLow(tLow>length(tInp))=length(tInp);

tHigh(tHigh==0) = 1;
tHigh(tHigh>length(tInp))=length(tInp);

for fNo = 1:1:length(fInp)
    %Now find the peaks
    SAIFFT(:,fNo) = SAIFFT(:,fNo).*winVal;
    [peakVal,peakLoc] = findpeaks(SAIFFT(:,fNo));
    %now take the biggest peak;
    [maxPeak,maxPeakInd] = max(peakVal);
    if(~isempty(maxPeak))
        peakInd = peakLoc(maxPeakInd(1,1),1);
        %peakInd = peakInd + tLow(fNo,1)-1;
    end
    
    if(~isempty(maxPeak))
        timePick(fNo,1) = tInp(1,peakInd);
        
        pickMag(fNo,1) = maxPeak(1,1);
        velPick(fNo,1) = stnDist/timePick(fNo,1);
        
    else
        timePick(fNo,1) = NaN;
        pickMag(fNo,1) = NaN;
        velPick(fNo,1) = NaN;
    end
end

velNew = velPick(~isnan(velPick),1);
fExtNew = fInp(~isnan(velPick),1);
% interpolate it on fExt
if(~isempty(velNew))
    velInterp = interp1(fExtNew,velNew,fInp);
    velPickUpdated = velInterp;
    timePickUpdated = stnDist./velPickUpdated;
else
    velPickUpdated = NaN*ones(length(fInp),1);
    timePickUpdated = NaN*ones(length(fInp),1);
end

end