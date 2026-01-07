function [SAIFFT,fExt,distOut,timePickUpdated,velPickUpdated,timePickResolve] = ...
                                          FTAN(inputCC,lagV,alpha,fSamp,fExt,n1,n2,plotSet)
%baser on Levshin & Ritzwoller 2001
%first perform fft of the data
[distOut] = calculateNodeDist(n1,n2);
tV = lagV/fSamp;

causalCC = inputCC(tV>=0,1);
acausalCC = inputCC(tV<0,1);
% 
if(max(abs(causalCC))>max(abs(acausalCC)))
    symmCC = causalCC;
else
    symmCC = [causalCC(1,1);flipud(acausalCC)];
end
% symmCC = (causalCC(2:end,1)+ flipud(acausalCC))/2;
% symmCC = [causalCC(1,1);symmCC];
tSymm = (1:1:length(symmCC))/fSamp;
%now get signal to noise ratio
[ccSNR,tSNR] = getSNR(symmCC,fSamp,5,9,0.2);
lSymmCC = length(symmCC);
%window the symmetric signal
symmCC = symmCC.*tukeywin(lSymmCC,0.01);
[fftCC] = fft(symmCC,lSymmCC);
fftCC = fftshift(fftCC);
fCC = fSamp/2*linspace(-1,1,lSymmCC);
% 
% subplot(3,2,3)
% plot(fCC,abs(fftCC));

SA = (1+sign(fCC))'.*fftCC;

tIFFT = (1:1:lSymmCC)/fSamp;
%now filter the amplitude spectra and transform it back to time
for fNo = 1:1:length(fExt)
    [filterVal] = getGaussFilter(fExt(fNo,1),alpha,fCC);
    SAFilt = SA.*filterVal';
    SAFilt = ifftshift(SAFilt);
    fCCNew = fSamp*linspace(0,1,lSymmCC);
    SAIFFT(:,fNo) = abs((ifft(SAFilt)));
    %Now find the peaks
    [peakVal,peakLoc] = findpeaks(SAIFFT(:,fNo));
    %now take the biggest peak;
    [maxPeak,maxPeakInd] = max(peakVal);
    peakInd = peakLoc(maxPeakInd,1);
    peakResolve = maxPeak*0.98;
    
    if(~isempty(maxPeak))
        timePick(fNo,1) = tIFFT(1,peakLoc(maxPeakInd(1,1),1));
        peakResolveInd =find(SAIFFT(peakInd:end,fNo)<peakResolve,1,'first');
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
        pickMag(fNo,1) = NaN;
        velPick(fNo,1) = distOut/timePick(fNo,1);
    end
end
%normalize the magnitude to 1
pickMag = pickMag./max(max(SAIFFT));
SAIFFTAbs = SAIFFT;
SAIFFT = SAIFFT/max(max(SAIFFT));

%now get rid of bad peaks
badVelInd = find(velPick>600 | velPick<50);
velPickUpdated = velPick;
velPickUpdated(badVelInd,1) = NaN;
%update again based on SNR
lowMagInd = find(pickMag<0.5);
velPickUpdated(lowMagInd,1) = NaN;

%now get rid of bad timePick
timePickUpdated = distOut./velPickUpdated;
timePickResolve(badVelInd,1) = NaN;
timePickResolve(lowMagInd,1) = NaN;
% subplot(3,2,5)
% plot(fCCNew,abs(SAFilt)); 
% subplot(3,2,6)
% plot(tIFFT,abs(SAIFFT));

if(plotSet == 1)
    figure(1)
    subplot(2,3,1)
    plot(tV,inputCC);
    xlim([-5,5]);
    hold on;
    plot(tSymm,symmCC,'r');
    hold off
    
    subplot(2,3,2)
    plot(fCC,abs(SA));
    
    subplot(2,3,3)
    imagesc(fExt,tIFFT,(SAIFFTAbs));
    shading interp
    colorbar;
    colormap('jet');
    hold on;
    plot3(fExt,timePickUpdated,pickMag,'k','LineWidth',2);
    plot3(fExt,timePickResolve,pickMag,'k','LineWidth',2);
    ylim([0,5]);
    hold off;
    
    subplot(2,3,4)
    plot(tSNR,ccSNR,'r','LineWidth',2);
    
    subplot(2,3,5)
    plot(fExt,velPickUpdated);
    hold on;
    plot(fExt,velPickUpdated,'ko','MarkerSize',6,'MarkerFaceColor','k');
    xlim([2,8]);
    ylim([0,800]);
    hold off;
%     
%     figure(100)
%     plotseis(SAIFFT,tIFFT,fExt,1);
end
end

