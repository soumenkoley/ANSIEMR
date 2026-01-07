% doFTAN
%clear;
%% run this code after assembleFiniteCCPairs.m
close all;
% this function performs FTAN on a desired station pair
%addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive2\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');
% load the picked phasse velocity for the subArray
fPathPhase = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\subArray00\';
load([fPathPhase,'phaseVelNewTaper.mat']);
% this file is from the small passive survey at Terziet
load('thGrpVelSmooth.mat');
load('HdBp3_8Hz.mat'); % a 3-8 Hz bandpass filter

%% parameters for FTAN
alphaVal = 300;
minAlpha = 30; maxAlpha = 50;
Beta = 3;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
n1 = "YHLWA"; n2= "Z2NAA";
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.1:0.05:8)';
plotSet = 1;
%% main code starts
% now, from the stnList, find the index of the stations
stnAInd = find(stnList==n1); stnBInd = find(stnList==n2);
n1Ind = find((rayAttribute(:,7)==stnAInd) | (rayAttribute(:,7)==stnBInd));
n2Ind = find((rayAttribute(:,8)==stnAInd) | (rayAttribute(:,8)==stnBInd));

pairInd = intersect(n1Ind,n2Ind);

if(~isempty(pairInd))
    nA = rayAttribute(pairInd,7);
    nB = rayAttribute(pairInd,8);
    
    if(nA == stnAInd)
        ccPair = ccStoreFinal(:,pairInd);
    else
        ccPair = flipud(ccStoreFinal(:,pairInd));
    end
else
    error('Station pair not found!');
end
% station pair read
%% perform FTAN
figure(90)
ccSee = ccStoreFinal(:,pairInd);
%ccSee = filtfilt(HdBp3_8Hz.Numerator,1,ccStoreFinal(:,pairInd));
plot(tArray,ccSee);
xlim([-10,10]);
lagVal = -500:1:500; % + 20 s and -20 s on ach side
tVec = tArray;
distNow = (rayAttribute(pairInd,5));

vMin = 0.5; vMax = 3; vInterval = 0.1;
%[w] = vel_taper2(tVec,1/fSamp,distNow/1000,vMin, vMax, vInterval );
%ccPair = ccPair.*(w');
%[newCC,tt] = doVelTaper(ccPair,fSamp,tVec,vMin,vMax,vInterval,distNow);
newCC = ccPair;
newCC = ccStoreFinal(:,pairInd);
%newCC = [0;diff(newCC)];

% load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\CCParamTest\NoInst\winSeg\ZRLWA-ZOKUA60min.mat');
% newCC = ccStore';
% newCC = newCC./max(abs(newCC));
%newCC = [0;diff(newCC)];

%ccPair = meanCC;
slopeNow = rayAttribute(pairInd,6);
disp(['Slope now = ', num2str(slopeNow)]);
caseString = ["causal";"acausal";"symm"];

for signalType = 1:1:length(caseString)
    caseType = caseString(signalType,1);
    for filterNo = 1:1:length(alphaVal)
        plotNum = (signalType-1)*3+1;
        
        [fExt,distOut,tPU,velPickUpdated,phPick,figOut] = ...
                  FTANNew(newCC,lagVal,minAlpha,maxAlpha,fSamp,fExt,distNow,...
                   plotSet,caseType,plotNum,Beta);
        %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
        disp(['Distance is = ', num2str(distOut)]);
        
        figure(104)
        hold on;
        plot(fExt,velPickUpdated,'-o','color',colorsAll(signalType,:));
        %plot(fExt,velPickUpdated,'color',colorsAll(signalType,:));
        plot(2.6:0.25:8,vgBC,'cyan','LineWidth',2)
        title(['Case = ', caseString(signalType,1)]);
        %figure(104);hold on;
        ylim([100,500]);
        xlim([2,8]);
    end
end

%% now process the group vel
% now get the phase velocity
% By intergration
deltaOmega = 0.2;
[f1,t1] = getpts(figOut);
velGrp = distOut./t1;
% now call the function to convert group velocity to phase velocity
[phaseVelOut1,fPhaseOut] = calculatePhaseVel(phPick,fExt,velGrp,f1,distOut);

figHand = openfig([fPathPhase,'subArray00NewCC.fig']);
hold on;
%plot3(fExt,distOut./(phPick./(2*pi*fExt)),1000*ones(length(fExt),1),'y--','LineWidth',2);
for k = 1:1:length(phaseVelOut1(1,:))
    plot3(fPhaseOut,phaseVelOut1(:,k),1000*ones(length(fPhaseOut),1),'LineWidth',2);
end
plot3(xi,yi,1000*ones(length(xi),1),'LineWidth',2);
plot3(f1,velGrp,1000*ones(length(f1),1),'k--','LineWidth',2);
ylim([0,5000]);
% convert obs phase velocity to group velocity
phaseInp = [xi,yi];
[vgOut] = getObsGrp(phaseInp);
plot3(vgOut(:,1),vgOut(:,2),1000*ones(length(vgOut(:,1)),1),'k','LineWidth',2);
%% now time to input the transition frequency, that is the frequency 
% from which the group velocity will be converted to phase
prompt = "Enter transition frequency?";
transFreq = input(prompt);
transFreq = [1.5;2;2.5;2.8];
for transNo = 1:1:length(transFreq)
    transFreqInd = find(xi>=transFreq(transNo,1),1,'first');
    vInit = yi(transFreqInd,1); fInit = xi(transFreqInd,1);
    % now interpolate the velGrp picked starting from transFreq to 8 Hz
    newFreq = transFreq(transNo):deltaOmega:max(f1);
    velGrpIntp = interp1(f1,velGrp,newFreq);
    %grpInv = cumsum(1./velGrpIntp);
    %wNum = 2*pi*fInit/vInit + 2*pi*deltaOmega*grpInv;
    grpInv = cumtrapz(1./velGrpIntp);
    wNum = 2*pi*fInit/vInit + 2*pi*deltaOmega*grpInv;
    phaseVelOut = 2*pi*newFreq./wNum;
    
    % now plot in on the same surface plot
    hold on;
    plot3(newFreq,phaseVelOut,1000*ones(length(newFreq),1),'m--','LineWidth',2);
    hold off;
end
xlim([1,7])
% N = 4;
% lambdaCorr = 0;
% phaseVel = (distOut*instOmega)./(phasePick+(instOmega.*tPU)-pi/4-(2*pi*N)-lambdaCorr);
%
% figure(2)
% plot(instOmega/2/pi,phaseVel);

% get the theoretical travel times
thTime = distOut./vgBC;
thFreq = 2.6:0.25:8;

figure(1)
subplot(2,3,4)
hold on;
plot3(thFreq,thTime,10*ones(length(thFreq),1),'LineWidth',4,'color','green');
hold off;
