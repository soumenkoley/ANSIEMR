% run FTAN
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
minAlpha = 300; maxAlpha = 500;
Beta = 3;
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.1:0.1:8)';
plotSet = 1;
%% main code starts

% station pair read
%% perform FTAN
lagVal = -500:1:500; % + 20 s and -20 s on ach side
tVec = tArray;
pairNo = 158;
distNow = distBin(1,pairNo);

vMin = 0.1; vMax = 2; vInterval = 0.05;
%[w] = vel_taper2(tVec,1/fSamp,distNow/1000,vMin, vMax, vInterval );
%ccPair = ccPair.*(w');
%load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\AzAvgCC\cc110-38.mat');
%ccPair = meanCC;
caseString = ["symm"];
ccPair = stackedTrace(:,pairNo);
figure(90)
%ccSee = filtfilt(HdBp3_8Hz.Numerator,1,ccStoreFinal(:,pairInd));
plot(tArray,ccPair);
xlim([-10,10]);
plotNum = 3+1;

[fExt,distOut,tPU,velPickUpdated,phPick,figOut] = ...
    FTANNew(ccPair,lagVal,minAlpha,maxAlpha,fSamp,fExt,distNow,...
    plotSet,caseString,plotNum,Beta);
