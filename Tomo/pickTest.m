% this script was written to pick the travel times between station pairs by
% using the theoretical travel times as a reference
clear; close all;

% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-corr and the theoretical travel times
load('grpTheoTimes.mat');

stnA = "Z2MDA"; stnB = "0GK6A";

% FTAN parameters
minAlpha = 30; maxAlpha = 50;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1:0.05:3)';
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

load('test.mat');
ccStoreFull = ccStoreFull';
stnDist  = 2720;
figure(316)
for s = 1:1:length(caseString)
    caseType = caseString(s,1);
    
    plotNum = (s-1)*3+1;
    
    [fExt,distOut,SAOut,tOut,symmCC,figOut] = ...
        FTANNewUse(ccStoreFull,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
    %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
    disp(['Distance is = ', num2str(distOut)]);
    
    
    % now get the SNR corresponding to these time-picks
   
    subplot(1,3,s)
    surf(fExt,tOut,SAOut);
    shading interp;
    colorbar; colormap('jet');
    ylim([0,8]);
    view(2);
    hold on;
    %plot3(fExt,tPick(:,s),1000*ones(length(fExt),1),'ko','MarkerFaceColor','k');
    %plot3(1.3:0.1:2.5,thTimeNow,1000*ones(length(thTimeNow),1),'m-','LineWidth',3);
    %hold off;
    
end
