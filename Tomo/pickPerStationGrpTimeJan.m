% this script was written to pick the travel times between station pairs by
% by comparing the travel time between the causal and the acausal
% cross-correlations; We allow a maximum of 20 % deviation between travel
% times; we will further test other percentages
clear; close all;

%% load input files
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');
% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-correlations and station names
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCFundUpdApril.mat');
% also load the vgBC from Terziet circular array campaign
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\vgBCTerziet.mat');
% load the theoretical travel times
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTimeFundApril.mat');
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelOvertLimits.mat');
load('HdBp1_5to5Hz.mat');
% load the velocity search band
load('vBand.mat');

stnA = "X0K9A"; stnB = "YHLWA";
% load the star velocity
%phStar = load('A:\TestInver\StarsApril\Star28\PhUpdAprilTheo.txt');
vCCBF = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray03vFK.mat');

starGrpVel = load('A:\TestInver\StarsApril\Star40\FundUpdApril.txt');
%% FTAN parametersP% FTAN parameters
minAlpha = 30; maxAlpha = 50;
%minAlpha = 1000; maxAlpha = 1200;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.1:0.05:5)';
fShow = 1.6:0.05:5;
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

%% locating the station pair from among all the correlations
% find the station indices and hence the correct correlation
stnAInd = find(allStn==stnA);
stnBInd = find(allStn==stnB);

stnANumInd = find(rayAttribute(:,7)==stnAInd);
stnBNumInd = find(rayAttribute(:,8)==stnBInd);

indAB = intersect(stnANumInd,stnBNumInd,'stable');
if(isempty(indAB))
    stnANumInd = find(rayAttribute(:,7)==stnBInd);
    stnBNumInd = find(rayAttribute(:,8)==stnAInd);
    indAB = intersect(stnANumInd,stnBNumInd,'stable');
    ccNow = flipud(ccStoreFinal(:,indAB));
    stnDist = rayAttribute(indAB,5);
    thTimeNow = thTime(:,indAB);
else
    ccNow = ccStoreFinal(:,indAB);
    stnDist = rayAttribute(indAB,5);
    thTimeNow = thTime(:,indAB);
end

% apply the filter
%ccNow = filtfilt(HdBp1_5to5Hz.Numerator,1,ccNow);
% updating this on Feb 05, 2025
%thTimeNow = stnDist./thGrpVel;
%% travel time based on the Terziet group velocity measurement

grpTimeIni(:,1) = vgBC(:,1); % frequecy Hz
grpTimeIni(:,2) = stnDist./vgBC(:,2); % travel time in secs

vDiff = 0.2; % allowing for a 20% devviation of the causal or acausal travel
% time to the symmetric travel times

%% perform FTAN and time picking

% peform the analysis for the symmetric correlation
caseType = caseString(3,1);

plotNum = 1;
    
[fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);

disp(['Distance is = ', num2str(distOut)]);

%% updated on July 14, 2025
fftCC = fft(symmCC);
fVec = linspace(0,1,501)*25;
J0 = besselj(0,(2*pi*vCCBF.vAll(:,1)*distOut)./(vCCBF.vAll(:,2)));

figure(1007);
plot(fVec,real(fftCC),'b');
hold on;
plot(vCCBF.vAll(:,1),J0,'r');
hold off;
xlim([1,5]);

vv = 2*pi*distOut*fVec'./(unwrap(angle(fftCC)));
figure(1008);
hold on;
plot(fVec,abs(vv))
xlim([1.5,5]); ylim([200,4000]);
% convert the star group vel to phase velocity
c0 = 4000;
fSt = 1.0;
k0 = 2*pi*fSt/c0;
fStInd = find(starGrpVel(:,1)==fSt,1,'first');
fAll = starGrpVel(fStInd:end,1);
k = 2*pi*cumtrapz(fAll(1:end),1./(1*starGrpVel(fStInd:end,2)));
kC = k0-k(1,1);
k = kC+k;
phVelOut = 2*pi*fAll(1:end)./k;
plot(fAll,phVelOut,'m');
plot(vCCBF.vAll(:,1),vCCBF.vAll(:,2),'k');
legend({'Now','fromStar','CCBF'});
hold off;

% now the function to identify the arrival times
%[tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);

[vPickSymm,tPickSymm] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
%[vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
% now get the SNR corresponding to these time-picks
spectralSNR(:,1) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymm(:,1));

figure(1)
subplot(1,3,1);
surf(fExt,tOut,SAOut); shading interp;
hold on;
view(2);colormap('jet');
plot3(fExt,tPickSymm,ones(length(fExt),1),'k','LineWidth',2);
plot3(starGrpVel(:,1),distOut./starGrpVel(:,2),ones(length(starGrpVel(:,1)),1),'m','LineWidth',2);
plot3(starGrpVel(:,1),distOut./(1.2*starGrpVel(:,2)),ones(length(starGrpVel(:,1)),1),'m--','LineWidth',2);
plot3(starGrpVel(:,1),distOut./(0.8*starGrpVel(:,2)),ones(length(starGrpVel(:,1)),1),'m--','LineWidth',2);

ylim([0,10]);
hold off;

% perform FTAN of the causal part of the cross-correlation
caseType = caseString(1,1);

plotNum = 2;
    
[fExt,distOut,SAOut,tOut,cauCC,figOutCau] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);

% pick the travel times that are close to the time picks for the symmetric
% correlations

[vPickCau,tPickCau] = pickTTCauAcauJan(SAOut,fExt,tOut,distOut,vPickSymm,vDiff);
% now get the SNR corresponding to these time-picks
spectralSNR(:,2) = getSpectSNR(cauCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickCau(:,1));

figure(1)
subplot(1,3,2);
surf(fExt,tOut,SAOut); shading interp;
hold on;
view(2);colormap('jet');
plot3(fExt,tPickCau,ones(length(fExt),1),'ko','LineWidth',2);
%plot3(starGrpVel(:,1),distOut./starGrpVel(:,2),ones(length(starGrpVel(:,1)),1),'m','LineWidth',2);

ylim([0,10]);
hold off;

% perform FTAN of the acausal part of the cross-correlation
caseType = caseString(2,1);

plotNum = 3;
    
[fExt,distOut,SAOut,tOut,AcauCC,figOutACau] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);

% pick the travel times that are close to the time picks for the symmetric
% correlations

[vPickAcau,tPickAcau] = pickTTCauAcauJan(SAOut,fExt,tOut,distOut,vPickSymm,vDiff);
% now get the SNR corresponding to these time-picks
spectralSNR(:,3) = getSpectSNR(AcauCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickAcau(:,1));

% perform FTAN again and try to get the phase velocity
% [fExt,distOut,SAOut,tOut,symmCC,phVel,figOutSymm] = ...
%         FTANNewPhase(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
%         plotSet,caseType,plotNum,Beta,tPickSymm);
%     
figure(1)
subplot(1,3,3);
surf(fExt,tOut,SAOut); shading interp;
hold on;
view(2);colormap('jet');
plot3(fExt,tPickAcau,ones(length(fExt),1),'ko','LineWidth',2);
%plot3(starGrpVel(:,1),distOut./starGrpVel(:,2),ones(length(starGrpVel(:,1)),1),'m','LineWidth',2);

ylim([0,10]);
hold off;

vAll = [vPickSymm,vPickCau,vPickAcau];
for i = 1:1:length(fExt)
    nanInd = isnan(sum(vAll(i,1)));
    if(~nanInd)
        vFinal(i,1) = fExt(i,1);
        vFinal(i,2) = vAll(i,1); % fill in with the symmetric velocity
    else
        vFinal(i,1) = fExt(i,1);
        vFinal(i,2) = NaN; % fill in with the symmetric velocity
    end
end

% now check if the these SNRs are less than 5

for i = 1:1:length(fExt)
    minSNR = min(spectralSNR(i,1));
    if(minSNR<5)
        vFinal(i,2) = NaN;
    end
end

figure(2);
hold on;

plot(fExt,vAll(:,1),'g','LineWidth',2);
plot(fExt,vAll(:,2),'b','LineWidth',2);
plot(fExt,vAll(:,3),'r','LineWidth',2);
plot(fExt,vAll(:,1),'go');
plot(fExt,vAll(:,2),'bo');
plot(fExt,vAll(:,3),'ro');
plot(fExt,vFinal(:,2),'k*');
plot(starGrpVel(:,1),starGrpVel(:,2),'m','LineWidth',2);
plot(fAll,stnDist./thTimeNow,'c','LineWidth',2);
plot(fAll,1.4*stnDist./thTimeNow,'c--','LineWidth',2);
plot(fAll,0.6*stnDist./thTimeNow,'c--','LineWidth',2);

set(gca,'YScale','log');
xlim([1,5]);
hold off;
