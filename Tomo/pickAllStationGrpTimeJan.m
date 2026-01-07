% this script was written to pick the travel times between station pairs by
% by comparing the travel time between the causal and the acausal
% cross-correlations; We allow a maximum of 20 % deviation between travel
% times; Theoretical group travel times are also used as a reference
% we perform the pick for all stations
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

% load the velocity search band
load('vBand.mat');

stnA = "0PM7A"; stnB = "0FL9A";
%% FTAN parameters
% FTAN parameters
minAlpha = 30; maxAlpha = 50;
%minAlpha = 1000; maxAlpha = 1200;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.1:0.05:5)';
fShow = 1.0:0.05:5;
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];
ifPlot = 0;

%% locating the station pair from among all the correlations
% find the station indices and hence the correct correlation
% stnAInd = find(allStn==stnA);
% stnBInd = find(allStn==stnB);
%
% stnANumInd = find(rayAttribute(:,7)==stnAInd);
% stnBNumInd = find(rayAttribute(:,8)==stnBInd);
%
% indAB = intersect(stnANumInd,stnBNumInd,'stable');
% if(isempty(indAB))
%     stnANumInd = find(rayAttribute(:,7)==stnBInd);
%     stnBNumInd = find(rayAttribute(:,8)==stnAInd);
%     indAB = intersect(stnANumInd,stnBNumInd,'stable');
%     ccNow = flipud(ccStoreFinal(:,indAB));
%     stnDist = rayAttribute(indAB,5);
%     thTimeNow = thTime(:,indAB);
% else
%     ccNow = ccStoreFinal(:,indAB);
%     stnDist = rayAttribute(indAB,5);
%     thTimeNow = thTime(:,indAB);
% end

nCC = length(ccStoreFinal(1,:));
%% travel time based on the Terziet group velocity measurement
%
% grpTimeIni(:,1) = vgBC(:,1); % frequecy Hz
% grpTimeIni(:,2) = stnDist./vgBC(:,2); % travel time in secs

vDiff = 0.2; % allowing for a 20% devviation of the causal or acausal travel
% time to the symmetric travel times

%% perform FTAN and time picking

for ccNo = 1:1:nCC
    ccNow = ccStoreFinal(:,ccNo);
    stnDist = rayAttribute(ccNo,5);
    % peform the analysis for the symmetric correlation
    caseType = caseString(3,1);
    
    plotNum = 1;
    %if(stnDist>2000)
    [fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
    
    disp(['Distance is = ', num2str(distOut)]);
    
    % now the function to identify the arrival times
    %[tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);
    thTimeNow = thTime(:,ccNo);
    thVel = stnDist./thTimeNow;
%     
    %% Comment here
    % updating this on Feb 05, 2025
    %thTimeNow = stnDist./thGrpVel;
    [vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
    
    % now get the SNR corresponding to these time-picks
    spectralSNR(:,1) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymm(:,1));
    if(ifPlot)
        figure(1)
        subplot(1,3,1);
        hold on;
        surf(fExt,tOut,SAOut); shading interp;
        view(2);colormap('jet');
        plot3(fExt,tPickSymm,ones(length(fExt),1),'k','LineWidth',2);
        ylim([0,10]);
        hold off;
    end
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
    
    if(ifPlot)
        figure(1)
        subplot(1,3,2);
        hold on;
        surf(fExt,tOut,SAOut); shading interp;
        view(2);colormap('jet');
        plot3(fExt,tPickCau,ones(length(fExt),1),'ko','LineWidth',2);
        ylim([0,10]);
        hold off;
    end
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
    
    if(ifPlot)
        figure(1)
        subplot(1,3,3);
        hold on;
        surf(fExt,tOut,SAOut); shading interp;
        view(2);colormap('jet');
        plot3(fExt,tPickAcau,ones(length(fExt),1),'ko','LineWidth',2);
        ylim([0,10]);
        hold off;
    end
    vAll = [vPickSymm,vPickCau,vPickAcau];
    vFinal = [];
    for i = 1:1:length(fExt)
        % change vAll(i,1) to vAll(i,:) to use symm times
        nanInd = isnan(sum(vAll(i,1)));
        if(~nanInd)
            vFinal(i,1) = fExt(i,1);
            vFinal(i,2) = vAll(i,1); % fill in with the symmetric velocity
            grpTimeStore(i,ccNo) = stnDist/vAll(i,1);
            tErr = abs(stnDist/vAll(i,2)-stnDist/vAll(i,3));
            %grpTimeErrStore(i,ccNo) = tErr;
            grpTimeErrStore(i,ccNo) = grpTimeStore(i,ccNo)/10;
            grpSNRStore(i,ccNo) = spectralSNR(i,1);
        else
            vFinal(i,1) = fExt(i,1);
            vFinal(i,2) = NaN; % fill in with the symmetric velocity
            grpTimeStore(i,ccNo) = NaN;
            grpTimeErrStore(i,ccNo) = NaN;
            grpSNRStore(i,ccNo) = NaN;
        end
    end
    
    % now check if these SNRs are less than 5
    
    for i = 1:1:length(fExt)
        minSNR = min(spectralSNR(i,1));
        if(minSNR<10)
            vFinal(i,2) = NaN;
        end
    end
    
    if(ifPlot)
        figure(2);
        hold on;
        
        plot(fExt,vAll(:,1),'g','LineWidth',2);
        plot(fExt,vAll(:,2),'b','LineWidth',2);
        plot(fExt,vAll(:,3),'r','LineWidth',2);
        plot(fExt,vAll(:,1),'go');
        plot(fExt,vAll(:,2),'bo');
        plot(fExt,vAll(:,3),'ro');
        plot(fExt,vFinal(:,2),'k*');
        %plot(fShow,stnDist./thTimeNow,'k','LineWidth',2);
        %plot(fShow,1.3*stnDist./thTimeNow(:,1),'k--','LineWidth',2);
        %plot(fShow,0.7*stnDist./thTimeNow(:,1),'k--','LineWidth',2);
        
        set(gca,'YScale','log');
        xlim([1,5]);
        hold off;
    end
    disp(['Station = ',num2str(ccNo),' done'])
    %end
    %% Uncomment here
    

end