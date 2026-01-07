% this script was written to pick the group velocities following FTAN clean
% and a picker algorithm is used
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

%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\vgBCTerziet.mat');
% load the theoretical travel times

fTh = 1.0:0.05:5;
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTimeFundApril.mat');
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelOvertLimits.mat');

% load the velocity search band
load('vBand.mat');
load('vBandJuly.mat');
load('vBandFundAug.mat');

%% FTAN parameters
minAlpha = 30; maxAlpha = 50;
%minAlpha = 1000; maxAlpha = 1200;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;

fExt = (0.1:0.05:5)';
%fExt = (logspace(log10(0.5),log10(5),20))';
fExtClean = (logspace(log10(0.5),log10(5),20))';
fShow = 1.6:0.05:5;
plotSet = 0;
ifPlot = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];
vBandUL(:,1) = interp1(vLow(:,1),vLow(:,2),fExtClean);
vBandUL(:,2) = interp1(vUp(:,1),vUp(:,2),fExtClean);
%% some picking parameters
opts = struct();
opts.vmin = vBand(:,3);
opts.vmax = vBand(:,2);
opts.smoothness_lambda = 1e7;
opts.post_sg_win = 11;
opts.post_sg_poly = 2;
opts.debug = true;

%% loop over all cross-correlations
lSymmCC = 501;
tIFFT = (0:1:(lSymmCC-1))/fSamp;
tV = (-500:500)/fSamp;
fCC = linspace(0,1,lSymmCC)*fSamp;
nF = length(fExtClean);
for fNo = 1:1:nF
    alpha(fNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(nF-1)*(fNo-1))-1;
end

alpha = fliplr(alpha);

for ccNo = 1:1:length(rayAttribute(:,1))
    disp(['Doing cc = ',num2str(ccNo)]);
    ccNow = ccStoreFinal(:,ccNo);
    stnDist = rayAttribute(ccNo,5);
    thTimeNow = thTime(:,ccNo);
    
    %% perform FTAN and time picking
    
    % peform the analysis for the symmetric correlation
    caseType = caseString(3,1);
    
    plotNum = 1;
    
    causalCC = ccNow(tV>=0,1);
    acausalCC = ccNow(tV<0,1);
    %now symmetrize the signal
    symmCC = (causalCC(2:end,1)+ flipud(acausalCC))/2;
    symmCC = [causalCC(1,1);symmCC];
    symmCCCopy = symmCC;
    fftCC = fft(symmCC);
    
    %%
    
    %now filter the amplitude spectra and transform it back to time
    
    ccClean = zeros(lSymmCC,1);
    
    for fNo = 1:1:length(fExtClean)
        [filterVal] = getGaussFilterFull(fExtClean(fNo,1),alpha(fNo),fCC,Beta);
        %plot(fCC,filterVal,'b');
        filterValFull = [filterVal(1:251),fliplr(filterVal(2:251))];
        filtFFT = fftCC.*filterValFull';
        ifftCC = ifft(filtFFT);
        [tW] = genWW(vBandUL(fNo,1),vBandUL(fNo,2),stnDist,fSamp,tIFFT);
        ccClean = ccClean + ifftCC.*(tW)';
    end
    
    %% perform FTAN of clean signal
    ccClean = [zeros(500,1);ccClean];
    [fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccClean,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
    [vPickSymmClean,tPickSymmClean] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
    
    out = pick_group_curve_dpUpd((SAOut(:,19:end))', fExt(19:end,1), tIFFT, distOut, opts);
    spectralSNRCleanN(:,1) = getSpectSNR(symmCCCopy,fExt,minAlpha,maxAlpha,Beta,fSamp,[tPickSymmClean(1:18,1);out.times_pick']);
    
    %%
    if(ifPlot)
        figure(3)
        hold on;
        surf(fExt,tOut,SAOut); shading interp;
        view(2);colormap('jet');
        plot3(fTh,thTimeNow,ones(length(fTh),1),'m--','LineWidth',2);
        plot3(fTh,1/0.7*thTimeNow,ones(length(fTh),1),'g--','LineWidth',2);
        plot3(fTh,1/1.3*thTimeNow,ones(length(fTh),1),'g--','LineWidth',2);
        plot3(fExt(19:end),out.times_pick,ones(length(fExt(19:end)),1),'r','LineWidth',2);
        ylim([0,10]);
        hold off;
    end
    %spectralSNRAll = [fExt,spectralSNRCleanN];
    grpTimeStore(:,ccNo) = out.times_pick';
    grpTimeErrStore(:,ccNo) = out.times_pick'/10;
    grpSNRStore(:,ccNo) = spectralSNRCleanN;
    
    %disp(['Distance is = ', num2str(stnDist)]);
    
end