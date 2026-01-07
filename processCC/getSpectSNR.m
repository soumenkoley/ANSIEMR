function [spectralSNR] = getSpectSNR(inputCC,fExt,minAlpha,maxAlpha,beta,fSamp,tPick)
    %this function computes te spectral SNR as a function of frequency of
    %the input cross correlation
    % a slightly updated version as compared to the function
    % getSpectralSNR.m
    % tPick is the time around which we get the max of the signal
    tNoiseStart = 15; tNoiseEnd = 18; % units in seconds
    lenF = length(fExt);
    spectralSNR = zeros(lenF,1);
    for bandNo = 1:1:lenF
        alpha(bandNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(lenF-1)*(bandNo-1))-1;
        [filtSig] = bandFilt(inputCC,fSamp,alpha(bandNo),beta,fExt(bandNo,1));
        
        if(~isnan(tPick(bandNo,1)) && tPick(bandNo)<20)
            [spectralSNR(bandNo,1)] = getRatio(filtSig,tNoiseStart,tNoiseEnd,fSamp,tPick(bandNo,1));
        else
            [spectralSNR(bandNo,1)] = NaN;
        end
    end
    
%     figure(110)
%     hold on;
%     plot(fExt,spectralSNR,'r*');
%     plot(fExt,spectralSNR,'r');
%     ylim([0,17]);
end

