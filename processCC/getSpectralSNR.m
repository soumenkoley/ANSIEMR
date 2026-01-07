function [spectralSNR] = getSpectralSNR(inputCC,fExt,minAlpha,maxAlpha,beta,fSamp,plotNumber)
    %this function computes te spectral SNR as a function of frequency of
    %the input cross correlation
    tNoiseStart = 10; tNoiseEnd = 15; % units in seconds
    lenF = length(fExt);
    spectralSNR = zeros(lenF,1);
    for bandNo = 1:1:lenF
        alpha(bandNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(lenF-1)*(bandNo-1))-1;
        [filtSig] = bandFilt(inputCC,fSamp,alpha(bandNo),beta,fExt(bandNo,1));
        [spectralSNR(bandNo,1)] = getRatio(filtSig,tNoiseStart,tNoiseEnd,fSamp);
    end
    
%     figure(110)
%     hold on;
%     plot(fExt,spectralSNR,'r*');
%     plot(fExt,spectralSNR,'r');
%     ylim([0,17]);
end

