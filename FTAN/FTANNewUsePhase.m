function [fExt,distOut,SAIFFTAbs,phaseSA,tIFFT,symmCC,figTT] = ...
    FTANNewUsePhase(inputCC,lagV,minAlpha,maxAlpha,fSamp,fExt,distIn,plotSet,caseNo,picNo,beta)
%based on Levshin & Ritzwoller 1989,2001
% updated April 26, 2022
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
%[spectralSNR] = getSpectralSNR(symmCC,fExt,minAlpha,maxAlpha,beta,fSamp,picNo);
lSymmCC = length(symmCC);
%window the symmetric signal

[fftCC] = fft(symmCC,lSymmCC);
fCC = fSamp*linspace(0,1,lSymmCC);
%

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
end

SAIFFTAbs = SAIFFT;

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
    view(2);
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
else
    figTT = 0;
end
end