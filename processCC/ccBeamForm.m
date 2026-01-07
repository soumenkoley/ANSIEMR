 % this script beamforms the cross correlations over subblocks of the entire
% array and gets the direction and slowness at a desired frequency

clear; close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');
% first load the cross-correlation pairs and the rayAttribute file
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray051112.mat';
%fPath = 'B:\LimburgBigSurvey1CC-Pair\SubArrayCC\subArray24.mat';
load(fPath);

fStart = 2.2:0.1:2.7; fEnd = fStart+0.1;
fSamp = 25;

%pArray = (1/1000):(1/80000):(1/400); %  for subgroup B1,B2 at low freqs upto 4 Hz
pArray = (1/5000):(1/80000):(1/400); % for subgrp B1,B2 at high freqs
%pArray = (1/450):(1/10000):(1/50); % for subgroup B2 % frwqs 5.5-7Hz
phiArray = 0:2:360;
vArray = 1./pArray;
% vMax = [[5000,5000,5000,5000,5000,5000];...
%     [2500,2500,2500,2500,2500,2500];...
%     [2000,2000,2000,2000,2000,2000];...
%     [1500,1500,1500,1500,1500,1500]];
% vMin = [[3000,3000,3000,3000,3000,3000];...
%     [2000,2000,2000,2000,2000,2000];...
%     [1600,1600,1600,1600,1600,1600];...
%     [600,600,600,600,600,600]];
vMax = [[4000,4000,4000,4000,4000,4000]];
vMin = [[1000,1000,1000,1000,1000,1000]];
for i = 1:1:length(fStart)
    for modeNo = 1:1:1
        vMaxInd(modeNo,i) = find(vArray<=vMax(modeNo,i),1,'first');
        vMinInd(modeNo,i) = find(vArray<=vMin(modeNo,i),1,'first');
    end
end

% for doughnut plot
pShow1 = (1/5000); pShow2 = 1/400;
p1Ind = find(pArray>=pShow1,1,'first');
p2Ind = find(pArray>=pShow2,1,'first');

pShow = pArray(p1Ind:p2Ind);

phiSum1 = 0; phiSum2 = 360;
phiSum1Ind = find(phiArray>=phiSum1,1,'first');
phiSum2Ind = find(phiArray>=phiSum2,1,'first');

% normalize the cross-correlations
for i = 1:1:length(ccStoreFinal(1,:))
    if(ccStoreFinal(1,i)~=0)
        ccStoreFinal(:,i) = ccStoreFinal(:,i)/max(abs(ccStoreFinal(:,i)));
    end
end

% % apply the velocity taper
% for i = 1:1:length(ccStoreFinal(1,:))
%      tapOut = vel_taper2(tArray,1/fSamp,rayAttribute(i,5)/1000,0.2,0.8,0.05);
%      ccStoreFinal(:,i) = ccStoreFinal(:,i).*tapOut';
% end

% now do the beamforming following ruigrok et al, 2017
% define the pArray and phiArray

% now do fft of the data
fftData = fft(ccStoreFinal,length(ccStoreFinal(:,1)),1);
lenData = length(ccStoreFinal(:,1));
fVec = fSamp*linspace(0,1,lenData);
pairDist = rayAttribute(:,5);
pairAz = rayAttribute(:,6);
for freqBand = 1:1:length(fStart)
    fStartInd = find(fVec>=fStart(1,freqBand),1,'first');
    fEndInd = find(fVec>=fEnd(1,freqBand),1,'first');
    
    fSumVec = fVec(fStartInd:fEndInd);
    
    bpStoreAvg = zeros(length(pArray),length(phiArray));
    bpStore = zeros(length(pArray),length(phiArray));
    
    for freqNo=1:1:length(fSumVec)
        fIntInd = fStartInd+(freqNo-1);
        XFData = fftData(fIntInd,:);
        fInt = fSumVec(1,freqNo);
        gg = sqrt(-1);
        
        for i = 1:1:length(phiArray)
            for j = 1:1:length(pArray)
                % now get the delay for all receiver pairs
                delayNow = (2*pi*fInt*pArray(1,j))*((pairDist').*(cosd((pairAz')-phiArray(1,i))));
                cmplxMult = exp(gg*delayNow);
                corrData = XFData.*cmplxMult;
                
                % now sum it and store the abs value
                sumCorrData = sum(corrData,2);
                bpStore(j,i) = abs(sumCorrData)/(length(corrData(1,:)));
            end
        end
        
        bpStoreAvg = bpStoreAvg+bpStore;
    end
    bpStoreAvg = bpStoreAvg/length(fSumVec);
    % make the bpStoreAvg into a long vector
    beamValCol = bpStoreAvg(:);
    
    [~] = plotDoughnut(bpStoreAvg(p1Ind:p2Ind,:),pShow,phiArray,fStart,freqBand);
    
    % sum bp across azimuth
    bpSum = mean(bpStoreAvg(:,phiSum1Ind:phiSum2Ind),2);
    
    % find the max bp between vMin and vMax
    vPick(freqBand,1) = mean(fSumVec);
    for modeNo = 1:1:1
        vSmall = vArray(vMaxInd(modeNo,freqBand):vMinInd(modeNo,freqBand));
        [maxVal,maxInd] = max(bpSum(vMaxInd(modeNo,freqBand):vMinInd(modeNo,freqBand),1));
        vPick(freqBand,modeNo+1) = vSmall(maxInd);
    end
    
    figure(3)
    hold on;
    %subplot(2,3,freqBand)
    plot(1./pArray,bpSum);
    
    figure(2)
    subplot(2,3,freqBand)
    imagesc(phiArray,pArray,bpStoreAvg);
    shading interp
    colorbar;
    %obsFInd = find(obsDC(:,1)>=fInt,1,'first');
    %obsSlowness = 1/obsDC(obsFInd,2);
    
    %obsSlowVec = obsSlowness*ones(length(phiArray),1);
    dummyVal = 10*ones(length(phiArray),1);
    hold on;
    %plot3(phiArray,obsSlowVec,dummyVal,'k--','LineWidth',2);
    hold off;
    title([num2str(fSumVec(1,1)),'-',num2str(fSumVec(1,end)),' Hz']);
end