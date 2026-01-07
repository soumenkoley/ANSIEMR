% this script was written to pick the travel times between station pairs by
% using the theoretical travel times as a reference
% this does it for all station pairs

clear; close all;

% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-corr and the theoretical travel times
load('grpTheoTimes.mat');

nCC = length(ccStoreFinal(1,:));
% FTAN parameters
minAlpha = 30; maxAlpha = 50;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1.25:0.05:3)';
plotSet = 1;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

fTh = 1.3:0.1:2.5;
for ccNo = 1:1:nCC
    
    ccNow = ccStoreFinal(:,ccNo);
    stnDist = rayAttribute(ccNo,5);
    
    % now proceed with FTAN of ccNow
    % first get the theoretical time
    thTimeNow = thTime(:,ccNo);
    % interpolate this time on fExt vector
    thTimeNowIntp = interp1(fTh,thTimeNow,fExt);
    
    tPick = [];tErrOut=[];
    spctralSNR = [];
    for s = 1:1:length(caseString)
        caseType = caseString(s,1);
        
        plotNum = (s-1)*3+1;
        
        [fExt,distOut,SAOut,tOut,symmCC,figOut] = ...
            FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
            plotSet,caseType,plotNum,Beta);
        %disp(['Distance is = ', num2str(distOut)]);
        
        % now the function to identify the arrival times
        [tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);
        
        % now get the SNR corresponding to these time-picks
        spectralSNR(:,s) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPick(:,s));
        
    end
    % now time to use the most closest travel time
    tFinalPick = [];
    tFinalError =[];
    finalSNR = [];
    for i = 1:1:length(fExt)
        % get the minimum error of the three cases, symm, acausal and causal
        b = (~isnan(tPick(i,:)));
        nanFlag  = sum(b);
        if(nanFlag>0)
            % means there is at least one not-Nan element
            [~,minInd] = min(tErrOut(i,:));
            tFinalPick(i,1) = tPick(i,minInd);
            %tFinalPick(i,1) = (tFinalPick(i,1)+thTimeNowIntp(i,1))/2;
            
            tFinalError(i,1) = abs(tFinalPick(i,1)-thTimeNowIntp(i,1));
            finalSNR(i,1) = spectralSNR(i,minInd);
        else
            tFinalPick(i,1) = NaN;
            tFinalError(i,1)= NaN;
            finalSNR(i,1) = NaN;
        end
    end
    grpTimeStore(:,ccNo) = tFinalPick;
    grpTimeErrStore(:,ccNo) = tFinalError;
    grpSNRStore(:,ccNo) = finalSNR;
    disp(['Station = ',num2str(ccNo)]);
end
% end of code