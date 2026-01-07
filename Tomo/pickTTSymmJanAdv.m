function [velPickUpdated,timePickUpdated] = pickTTSymmJanAdv(SAIFFT,fInp,fShow,tInp,stnDist,thTime,fSamp)
% this function applies findpeaks to find the peaks in the FTAN image
% give out the velocities for the symmetric cross-correlation

% generate a tukey window to taper the end of SAIFFT
tkySmall = tukeywin(100,1);
winVal = [ones(451,1);tkySmall(51:100,1)];
lentInp = length(tInp);

velApp = stnDist./thTime;
velHigh = 1.2*velApp(:,1);
velLow = 0.8*velApp(:,1);

thTimeUp = stnDist./velLow;
thTimeLow = stnDist./velHigh;

%thTimeUp = stnDist./100*ones(length(fShow),1);
%thTimeLow = stnDist./3000*ones(length(fShow),1);

fStartInd = find(fInp>=fShow(1),1,'first');
fEndInd = find(fInp>=fShow(end),1,'first');
for fNo = 1:1:length(fShow)
    %Now find the peaks
    thTimeUpInd = find(tInp<=thTimeUp(fNo,1),1,'last');
    thTimeLowInd = find(tInp<=thTimeLow(fNo,1),1,'last');
    
    if(isempty(thTimeLowInd) || thTimeLowInd==lentInp)
        % the faster velocity provides a travel time band beyond 20 s
        % not possible to have a pick
        timePick(fNo,1) = NaN;
        pickMag(fNo,1) = NaN;
        velPick(fNo,1) = NaN;
    else
        if(isempty(thTimeUpInd))
            thTimeUpInd = lentInp;
        end
        if((thTimeUpInd-thTimeLowInd)<3)
            thTimeLowInd = thTimeLowInd-3;
            if(thTimeLowInd<1)
                thTimeLowInd = 1;
                if((thTimeUpInd-thTimeLowInd)<3)
                    thTimeUpInd = thTimeUpInd+3;
                end
            end
        end
        SAIFFT(:,fStartInd+fNo-1) = SAIFFT(:,fStartInd+fNo-1).*winVal;
        [peakVal,peakLoc] = findpeaks(SAIFFT(thTimeLowInd:thTimeUpInd,fStartInd+fNo-1));
        % possible that there are multiple peaks in this band
        % use the peak with maximum amplitude

        [maxPeak,maxPeakInd] = max(peakVal);
        if(~isempty(maxPeak))
            peakInd = thTimeLowInd + peakLoc(maxPeakInd(1,1),1)-1;
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
    
end
velPick(velPick<velLow,1)= NaN;
velPick(velPick>velHigh,1) = NaN;

velNew = NaN*ones(length(fInp),1);
velNew(fStartInd:fEndInd,1) = velPick;
%velNew = velPick(~isnan(velPick),1);
%fExtNew = fShow(~isnan(velPick));
% interpolate it on fExt
%velInterp = interp1(fExtNew,velNew,fInp);
velPickUpdated = velNew;
timePickUpdated = stnDist./velPickUpdated;

end