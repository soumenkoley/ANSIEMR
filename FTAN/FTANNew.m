function [fExt,distOut,timePickUpdated,velPickUpdated,phasePick,figTT] = ...
    FTANNew(inputCC,lagV,minAlpha,maxAlpha,fSamp,fExt,distIn,plotSet,caseNo,picNo,beta)
%based on Levshin & Ritzwoller 1989,2001
%first perform fft of the data
load('thGrpVelSmooth.mat');

[distOut] = distIn;
tV = lagV/fSamp;

causalCC = inputCC(tV>=0,1);
acausalCC = inputCC(tV<0,1);
%now symmetrize the signal

switch caseNo
    case 'symm'
        symmCC = (causalCC(2:end,1)+ flipud(acausalCC))/2;
        symmCC = [causalCC(1,1);symmCC];
    case 'acausal'
        symmCC = [causalCC(1,1);flipud(acausalCC)];
    case 'causal'
        symmCC = causalCC;
end

tSymm = (1:1:length(symmCC))/fSamp;

%now get signal to noise ratio
tN1 = 10; tN2 = 20;
%[spectralSNR] = getSpectralSNR(symmCC,fExt,alpha,fSamp,picNo);
lSymmCC = length(symmCC);
%window the symmetric signal
% now do a velocity taper

%[velFilt] = vel_taper2( tV,dt,distOut,vMin,vMax,vInterval);
%velFilt = velFilt(tV>=0);
%symmCC = symmCC.*(velFilt');

[fftCC] = fft(symmCC,lSymmCC);
fCC = fSamp*linspace(0,1,lSymmCC);
%fftCC = -fftCC.*((sqrt(-1))*2*pi*fCC');
%
% subplot(3,2,3)
% plot(fCC,abs(fftCC));
% set the frequencies greater than nyquist freq, = 0
SA = (1+sign(fCC))'.*fftCC;
SA((fCC>fSamp/2),1) = 0;

tIFFT = (1:1:lSymmCC)/fSamp;
%now filter the amplitude spectra and transform it back to time
nF = length(fExt);

for fNo = 1:1:nF
    alpha(fNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(nF-1)*(fNo-1))-1;
end
alpha = fliplr(alpha);

for fNo = 1:1:length(fExt)
    [filterVal] = getGaussFilter(fExt(fNo,1),alpha(fNo),fCC,beta);
    SAFilt = SA.*filterVal';
    fCCNew = fSamp*linspace(0,1,lSymmCC);
    SAIFFT(:,fNo) = abs((ifft(SAFilt)));
    phaseSA(:,fNo) = (angle(ifft(SAFilt)));
    %Now find the peaks
    [peakVal,peakLoc] = findpeaks(SAIFFT(:,fNo));
    %now take the biggest peak;
    [maxPeak,maxPeakInd] = max(peakVal);
    if(~isempty(maxPeak))
        peakInd = peakLoc(maxPeakInd(1,1),1);
    end
    peakResolve = maxPeak*0.95;
    
    if(~isempty(maxPeak))
        timePick(fNo,1) = tIFFT(1,peakLoc(maxPeakInd(1,1),1));
        peakResolveInd =find(SAIFFT(peakInd:end,fNo)<peakResolve,1,'first');
        phasePick(fNo,1) = phaseSA(peakInd,fNo);
        %instOmega(fNo,1) = (phaseSA(peakInd+1,fNo)-phaseSA(peakInd,fNo))*fSamp;
        if(~isempty(peakResolveInd))
            timePickResolve(fNo,1) = tIFFT(1,(peakInd+peakResolveInd(1,1)-1));
        else
            timePickResolve(fNo,1) = NaN;
        end
        
        pickMag(fNo,1) = maxPeak(1,1);
        velPick(fNo,1) = distOut/timePick(fNo,1);
        %disp(['fNo = ', num2str(fNo)]);
        %SAIFFT(:,fNo) = SAIFFT(:,fNo)/max(SAIFFT(:,fNo));
    else
        timePick(fNo,1) = NaN;
        %phasePick(fNo,1) = NaN;
        pickMag(fNo,1) = NaN;
        velPick(fNo,1) = distOut/timePick(fNo,1);
        %instOmega(fNo,1) = NaN;
    end
end
%normalize the magnitude to 1
pickMag = pickMag./max(max(SAIFFT));
SAIFFTAbs = SAIFFT;
SAIFFT = SAIFFT/max(max(SAIFFT));

%now get rid of bad peaks
badVelInd = find(velPick>500 | velPick<50);
velPickUpdated = velPick;
%velPickUpdated(badVelInd,1) = NaN;
%update again based on SNR
lowMagInd = find(pickMag<0.3);
%velPickUpdated(lowMagInd,1) = NaN;

% now make a new array with non-nan values
velNew = velPickUpdated(~isnan(velPickUpdated),1);
fExtNew = fExt(~isnan(velPickUpdated),1);
% interpolate it on fExt
velInterp = interp1(fExtNew,velNew,fExt);
%velInterp = smooth(velInterp);
velPickUpdated = velInterp;
%now get rid of bad timePick
timePickUpdated = distOut./velPickUpdated;
%timePickResolve(badVelInd,1) = NaN;
%timePickResolve(lowMagInd,1) = NaN;

%errorTimePick = abs(timePickUpdated-timePickResolve);
%phasePick(badVelInd,1) = NaN;
%phasePick(lowMagInd,1) = NaN;

%instOmega(badVelInd,1) = NaN;
%instOmega(lowMagInd,1) = NaN;
% subplot(3,2,5)
% plot(fCCNew,abs(SAFilt));
% subplot(3,2,6)
% plot(tIFFT,abs(SAIFFT));

if(plotSet == 1)
    figure(picNo)
    subplot(2,3,1)
    plot(tV,inputCC);
    xlim([-10,10]);
    
    subplot(2,3,2)
    plot(tSymm,symmCC,'r');
    xlim([-10,10]);
    
    subplot(2,3,3)
    plot(fCC,abs(SA));
    
    subplot(2,3,4)
    surf(fExt,tIFFT,(SAIFFTAbs));
    shading interp
    colorbar;
    colormap('jet');
    hold on;
    plot3(fExt,timePickUpdated,pickMag,'k','LineWidth',2);
    plot3(fExt,timePickResolve,pickMag,'k--','LineWidth',2);
    ylim([0,9]);
    hold off;
    
    %     subplot(2,3,5)
    %     plot(tSNR,ccSNR,'r','LineWidth',2);
    
    subplot(2,3,6)
    plot(fExt,(velPickUpdated));
    hold on;
    plot(2.6:0.25:8,vgBC,'r','LineWidth',2)
    vgUp = 1.35*vgBC; vgLow = 0.65*vgBC;
    plot(2.6:0.25:8,vgUp,'r--','LineWidth',2);
    plot(2.6:0.25:8,vgLow,'r--','LineWidth',2);
    plot(fExt,velPickUpdated,'ko','MarkerSize',6,'MarkerFaceColor','k');
    xlim([2,8]);
    ylim([0,1000]);
    hold off;
    
    figure(2)
    imagesc(fExt,tIFFT,(SAIFFTAbs));
    shading interp
    colorbar;
    colormap('jet');
    hold on;
    plot3(fExt,timePickUpdated,pickMag,'k','LineWidth',2);
    plot3(fExt,timePickResolve,pickMag,'k--','LineWidth',2);
    ylim([0,6]);
    hold off;
    
    figTT = figure(66);
    surf(fExt,tIFFT,(SAIFFTAbs));
    shading interp
    colorbar;
    colormap('jet');
    hold on;
    plot3(fExt,timePickUpdated,pickMag,'k','LineWidth',2);
    plot3(fExt,timePickResolve,pickMag,'k--','LineWidth',2);
    ylim([0,9]);
    view(2);
    hold off;
    
    
end
end