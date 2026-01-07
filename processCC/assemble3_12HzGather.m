% this script was written to see the cross-correlation gather
% corresponding to the 3.125 Hz noise
% all cross-correlations are loaded and then we use only those
% that have a spectral SNR greater than 5
% azimuthal averaging is performed
% lets see how it goes

clear; close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');

% load station locations
load('nodeLocationsCartesian.mat');
%load all ccs
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCsNew.mat');
% load the filter
load('HdBp3_5Hz.mat');
fSamp = 25;
minDist = 25; maxDist = 8000;
fVec = linspace(0,1,1001)*fSamp;
f1Peak = 3.05; f2Peak = 3.2;
f1Floor = 3.2; f2Floor = 3.5;
f1PeakInd = find(fVec>=f1Peak,1,'first');
f2PeakInd = find(fVec>=f2Peak,1,'first');
f1FloorInd = find(fVec>=f1Floor,1,'first');
f2FloorInd = find(fVec>=f2Floor,1,'first');

% remove NaN cross-corrs
goodInd = (~isnan(ccStoreFinal(1,:)));
ccStoreFinal = ccStoreFinal(:,goodInd);
rayAttribute = rayAttribute(goodInd,:);
goodInd = [];
% now do one thing check the fft amplitudes and if they are all balanced
%fftCCAll = fft(ccStoreFinal,1001,1);

% have to loop over all traces for filtering
distBin = minDist:minDist:maxDist;
nTraces = length(ccStoreFinal(1,:));
avgGather = zeros(length(tArray),length(distBin));
gatherCount = zeros(length(distBin),1);
vel = zeros(nTraces,1);

figure(1);
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'b^');
hold on;

velMin = 1900; velMax = 3100;
colVec = jet(100);
for i = 1:1:nTraces
    fftNow = fft(ccStoreFinal(:,i));
    pkEnergy = max(abs(fftNow(f1PeakInd:f2PeakInd,1)));
    floorEnergy = mean(abs(fftNow(f1FloorInd:f2FloorInd,1)));
    snr(i,1) = pkEnergy/floorEnergy;
    if(snr(i,1)>=3)
        % good trace, filter it now
        ccFilt = filtfilt(HdBp3_5Hz.Numerator,1,ccStoreFinal(:,i));
        ccFilt = ccFilt/(max(abs(ccFilt)));
        distInd = floor(rayAttribute(i,5)/minDist) + 1;
        %pkTimeA = rayAttribute(i,5)/2000; pkTimeB = rayAttribute(i,5)/3000;
        [maxVal,maxInd] = max(abs(ccFilt));
        pkTime = tArray(maxInd);
        
        vel(i,1) = rayAttribute(i,5)/abs(pkTime);
        if((vel(i,1)>2000)&&(vel(i,1)<3000))
            avgGather(:,distInd) = avgGather(:,distInd) + ccFilt;
            gatherCount(distInd,1) = gatherCount(distInd,1) + 1;
            startXY = [rayAttribute(i,1),rayAttribute(i,3)];
            endXY = [rayAttribute(i,2),rayAttribute(i,4)];
            colInd = floor((vel(i,1)-velMin)/(velMax-velMin)*100);
            plot(startXY,endXY,'-','color',colVec(colInd,:));
        end
        
    end
end
hold off;
% normalize the average gather
for i = 1:1:length(distBin)
    if(gatherCount(i,1)~=0)
        avgGather(:,i) = avgGather(:,i)/gatherCount(i,1);
        avgGather(:,i) = avgGather(:,i)/max(abs(avgGather(:,i)));
    end
   
end
figure(30);
plotseis(avgGather,tArray,(distBin),1);

figure(31)
imagesc(distBin,tArray,avgGather)
colorbar;colormap('jet');
%camroll(-90);

% distBin = [distBin,max(rayAttribute(:,5))];
% 
% figure(1)
% 
% % now sum the traces in distance bins
% 
% nBins = length(distBin(1,:));
% for i = 1:1:(nBins-1)
%     [ind,~] = find(rayAttribute(:,5) >= distBin(i) & rayAttribute(:,5)<distBin(i+1));
%     stackedTrace(:,i) = mean(ccStoreFilt(:,ind),2);
%     stackedTrace(:,i) = stackedTrace(:,i)/max(abs(stackedTrace(:,i)));
% end
% figure(30);
% plotseis(stackedTrace,tArray,(distBin(1:(end-1))/1000),1);
% camroll(-90);
% 
% figure(31)
% imagesc((distBin(1:(end-1))/1000),tArray,stackedTrace);

