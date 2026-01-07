% this script was writen to see which correlations have velocity around 500
% m/s in the frequency band 2-5 Hz

clear; close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
% load all the correlations
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCs.mat');

% load the station locations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% load the filter 2_5 Hz
load('HdBp2_5Hz.mat');

nCC = length(rayAttribute(:,1));
tN1 = 15; tN2 = 20;
tN1Ind = find(tArray>=tN1,1,'first');
tN2Ind = find(tArray>=tN2,1,'first');

% loop across correlations
useCC = 1;
for i = 1:1:nCC
    filtCC = filtfilt(HdBp2_5Hz.Numerator,1,ccStoreFinal(:,i));
    rmsNoise = rms(filtCC(tN1Ind:tN2Ind,1));
    % get the hlbert transfor for signal envelope
    envSig = abs(hilbert(filtCC));
    [maxVal,maxInd] = max(envSig);
    snrPeak = maxVal(1,1)/rmsNoise;
    appVel = rayAttribute(i,5)/abs(tArray(maxInd));
    if(appVel<500 && appVel>100 && snrPeak>20)
        % get these correlations
        ccStore(:,useCC) = filtCC;
        rayStore(useCC,:) = rayAttribute(i,:);
        useCC = useCC+1;
    end
    %disp('one done');
end

% sort the CCs in 100 m bins
% now sum the traces in distance bins
distBin = min(rayStore(:,5)):50:max(rayStore(:,5));
distBin = [distBin,max(rayStore(:,5))];
nBins = length(distBin(1,:));

for i = 1:1:(nBins-1)
    [ind,~] = find(rayStore(:,5) >= distBin(i) & rayStore(:,5)<distBin(i+1));
    stackedTrace(:,i) = mean(ccStore(:,ind),2,'omitnan');
    stackedTrace(:,i) = stackedTrace(:,i)/max(abs(stackedTrace(:,i)));
end
figure(30);
plotseis(stackedTrace,tArray,(distBin(1:(end-1))/1000),1);
camroll(-90);

