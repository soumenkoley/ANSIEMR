function [snrOut] = getRatio(filtSignal,tNoiseStart,tNoiseEnd,fSamp,tPeak)
    %this function computes the ratio of the peak of the CC to the rms of
    %the noise band;
    dt = 1/fSamp;
    tSigStart = tPeak-0.2;
    tSigEnd = tPeak+0.2;
    
    if(tSigStart<=0)
        tSigStart = 0;
    end
    if(tSigEnd>20)
        tSigEnd=20;
    end
    sigInd = (tSigStart:dt:tSigEnd)*fSamp + 1;
    sigInd = round(sigInd);
    if(max(sigInd)>501)
        disp('Stop');
    end
    peakSig = max(abs(filtSignal(sigInd,1)));
    
    noiseInd = (tNoiseStart:dt:tNoiseEnd)*fSamp + 1;
    noiseSignal = filtSignal(round(noiseInd),1);
    rmsNoise = rms(noiseSignal);
    
    snrOut = peakSig/rmsNoise;
end

