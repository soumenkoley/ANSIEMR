% this script beamforms the cross correlations over subblocks of the entire
% array and gets the direction and slowness at a desired frequency

clear; %close all;

% first load the cross-correlation pairs and the rayAttribute file
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray0511.mat';

% load the vFK velocity for this subArray
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray05vFK.mat');

%fPath = 'B:\LimburgBigSurvey1CC-Pair\SubArrayCC\subArray04.mat';
load(fPath);
load('Fund.txt');
%ccStoreFinal = linearStack;
fStart = 1:0.1:3.5; fEnd = fStart;
fSamp = 25;

%pArray = (1/1000):(1/20000):(1/400); %  for subgroup B1,B2 at low freqs upto 4 Hz
pArray = (1/4000):(1/40000):(1/400); % for subgrp B1,B2 at high freqs
% pArray = (1/450):(1/10000):(1/100); % for subgroup B2 % frwqs 5.5-7Hz
phiArray = 0:1:359;
phi1 = 0; phi2 = 359;
phi1Ind = find(phiArray>=phi1,1,'first');
phi2Ind = find(phiArray>=phi2,1,'first');

% now do the beamforming following ruigrok et al, 2017
% define the pArray and phiArray
% load('HdBp1_5to4Hz.mat')
% xx = zeros(500,length(ccStoreFinal(1,:)));
% ccStoreFinal = [xx;ccStoreFinal];
% ccStoreFinal = [ccStoreFinal;xx];
% % now filter and see which correlations are symmetric
% for i = 1:1:length(ccStoreFinal(1,:))
%     ccStoreFinal(:,i) = filtfilt(HdBp1_5to4Hz.Numerator,1,ccStoreFinal(:,i));
% end
%tArray = (-1000:1:1000)/fSamp;
% normalize each trace

for i = 1:1:length(ccStoreFinal(1,:))
    if(ccStoreFinal(1,i)~=0)
        ccStoreFinal(:,i) = ccStoreFinal(:,i)/max(abs(ccStoreFinal(:,i)));
    end
end

% % apply the velocity taper
% for i = 1:1:length(ccStoreFinal(1,:))
%     tapOut = vel_taper2(tArray,1/fSamp,rayAttribute(i,5)/1000,0.15,1,0.05);
%     ccStoreFinal(:,i) = ccStoreFinal(:,i).*tapOut';
% end

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
    
    % sum the beampower across azimuth
    sumBP(:,freqBand) = mean(bpStoreAvg(:,phi1Ind:phi2Ind),2);
    sumBP(:,freqBand) = sumBP(:,freqBand)/max(sumBP(:,freqBand));
    
    disp(['Freq Band = ',num2str(fStart(freqBand))]);
    %figure(2)
    %subplot(2,3,freqBand)
    %imagesc(phiArray,pArray,bpStoreAvg);
    %shading interp
    %colorbar;
    %obsFInd = find(obsDC(:,1)>=fInt,1,'first');
    %obsSlowness = 1/obsDC(obsFInd,2);
    
    %obsSlowVec = obsSlowness*ones(length(phiArray),1);
    dummyVal = 10*ones(length(phiArray),1);
    %hold on;
    %plot3(phiArray,obsSlowVec,dummyVal,'k--','LineWidth',2);
    %hold off;
    %title([num2str(fSumVec(1,1)),'-',num2str(fSumVec(1,end)),' Hz']);
end

figHand = figure(78);
hold on;
surf((fStart+fEnd)/2,1./pArray,sumBP);
shading interp;
view(0,90);
%[fVal,phVel] = getpts(figHand);
phVelPath =  'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayPhVel\';
load([phVelPath,'SubArray2MelleschetFund.mat']);
plot3(fVal,phVel,1000*ones(length(fVal),1),'k');
load([phVelPath,'SubArray2MelleschetOvert1.mat']);
plot3(fVal,phVel,1000*ones(length(fVal),1),'k--');
%plot3(Fund(:,1),Fund(:,2),1000*ones(length(Fund(:,1)),1),'k','LineWidth',2);
xlim([1,8]);
plot3(vAll(:,1),vAll(:,2),1000*ones(length(vAll),1),'k');
hold off;
colorbar; colormap('jet');

% also now show the nodes in the subArray
A = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
refStn = "WUJ2A";
refLong = 5.9060966; refLat = 50.7513927;

refLatLong = [refLat,refLong];
R = 6371000; % in meters
for i = 1:1:length(allStn)
    endStnInd = find(allStn==allStn(i,1),1,'first');
    latLongNow = A(endStnInd,1:2);
    [xlong,ylat] = calculatedist(latLongNow,refLatLong,R);
    nodeLocCart(i,1:3) = [xlong,ylat,A(endStnInd,3)];
end

R = 6371000; % in meters
for i = 1:1:length(stnList)
    subArrayStnInd = find(allStn==stnList(i,1),1,'first');
    latLongNow = A(subArrayStnInd,1:2);
    [xlong,ylat] = calculatedist(latLongNow,refLatLong,R);
    nodeLocCartSubArray(i,1:2) = [xlong,ylat];
end

xAlt = min(nodeLocCart(:,1)):200:max(nodeLocCart(:,1));
yAlt = min(nodeLocCart(:,2)):200:max(nodeLocCart(:,2));
[xyAltA,xyAltB] = meshgrid(xAlt,yAlt);
% xyAltA = xyAltA(:);
% xyAltB = xyAltB(:);

% now get the altitude map
altIntp = griddata(nodeLocCart(:,1),nodeLocCart(:,2),nodeLocCart(:,3),...
    xyAltA,xyAltB);
altIntpMat = vec2mat(altIntp,length(xAlt));

figure(41)
h = imagesc(xAlt,yAlt,altIntpMat);
set(h, 'AlphaData', ~isnan(altIntpMat))
set(gca,'YDir','normal')
colorbar; colormap('jet');
hold on;
plot(nodeLocCart(:,1),nodeLocCart(:,2),'mo','MarkerSize',8,'MarkerFaceColor','k');
plot(nodeLocCartSubArray(:,1),nodeLocCartSubArray(:,2),'ko','MarkerSize',8,'MarkerFaceColor','r');
hold off;