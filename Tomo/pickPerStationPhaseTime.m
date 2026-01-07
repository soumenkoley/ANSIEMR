% this script was written to pick the travel times between station pairs by
% using the theoretical travel times as a reference
%clear; %close all;
clear;
fP = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\';
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');
% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-corr and the theoretical travel times
load('grpTheoTimes.mat');

stnA = "YYMLA"; stnB = "X4MBA";
% FTAN parameters
minAlpha = 30; maxAlpha = 50;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1.25:0.05:10)';
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

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
else
    ccNow = ccStoreFinal(:,indAB);
    stnDist = rayAttribute(indAB,5);
end
%ccNow = ccStoreFull';
%stnDist = 200;
%%ccFilt = filtfilt(F,1,ccNow);

% now proceed with FTAN of ccNow
% first get the theoretical time
thTimeNow = thTime(:,indAB);
fTh = 1.3:0.1:2.5;
% interpolate this time on fExt vector
thTimeNowIntp = interp1(fTh,thTimeNow,fExt);
figure(316)
for s = 1:1:length(caseString)
    caseType = caseString(s,1);
    
    plotNum = (s-1)*3+1;
    
    [fExt,distOut,SAOut,phOut,tOut,symmCC,figOut] = ...
        FTANNewUsePhase(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
    %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
    disp(['Distance is = ', num2str(distOut)]);
    
    % now the function to identify the arrival times
    %[tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);
    [tPick(:,s),phPick(:,s)] = pickTTPhase(SAOut,phOut,fExt,tOut);
    % now get the SNR corresponding to these time-picks
    
    %spectralSNR(:,s) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPick(:,s));
    figure(316)
    subplot(1,3,s)
    surf(fExt,tOut,SAOut);
    shading interp;
    colorbar; colormap('jet');
    ylim([0,8]);
    view(2);
    hold on;
    plot3(fExt,tPick(:,s),1000*ones(length(fExt),1),'ko','MarkerFaceColor','k');
    plot3(1.3:0.1:2.5,thTimeNow,1000*ones(length(thTimeNow),1),'m-','LineWidth',3);
    hold off;
    
    % now do the conversion from group to phase velocity
    grpVel = distOut./tPick(:,s);
    [phaseVelOut,fOut] = calculatePhaseVel(phPick(:,s),fExt,grpVel,fExt,distOut);
    disp('Trying')
    
    figure(21)
    subplot(1,3,s)
    
    hold on;
    plot(fOut,phaseVelOut);
    ylim([0,2000]);
    
    figure(22)
    subplot(1,3,s)
    hold on;
    plot(fExt,grpVel);
    ylim([0,2000])
end
% now time to use the most closest travel time
