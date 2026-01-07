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
tPos = tArray(tArray>=0);
tNeg = tArray(tArray<=0);

minVel = 100; maxVel = 600; %  units in m/s
ifFK = 0;
% loop across correlations
useCC = 1;
for i = 1:1:nCC
    filtCC = filtfilt(HdBp2_5Hz.Numerator,1,ccStoreFinal(:,i));
    rmsNoise = rms(filtCC(tN1Ind:tN2Ind,1));
    negCC = filtCC(tArray<=0,1);
    posCC = filtCC(tArray>=0,1);
    
    % get the hlbert transfor for signal envelope
    envSigNeg = abs(hilbert(negCC));
    [maxVal,maxInd] = max(envSigNeg);
    snrPeakNeg = maxVal(1,1)/rmsNoise;
    appVelNeg = rayAttribute(i,5)/abs(tNeg(maxInd));
    
    envSigPos = abs(hilbert(posCC));
    [maxVal,maxInd] = max(envSigPos);
    snrPeakPos = maxVal(1,1)/rmsNoise;
    appVelPos = rayAttribute(i,5)/abs(tPos(maxInd));
    
    % check how symmetric it is
    appVelDiff = abs(appVelPos-appVelNeg);
    appVelDiffPerct = appVelDiff/max(appVelPos,appVelNeg)*100;
    
    snrPeak = max(snrPeakPos,snrPeakNeg);
    
    appVelPosFlag = (appVelPos<maxVel && appVelPos>minVel);
    appVelNegFlag = (appVelNeg<maxVel && appVelNeg>minVel);
    
    if(appVelPosFlag && appVelNegFlag && appVelDiffPerct<20 && snrPeak>10)
        % get these correlations
        ccStore(:,useCC) = filtCC;
        rayStore(useCC,:) = rayAttribute(i,:);
        useCC = useCC+1;
    end
    %disp('one done');
end

% sort the CCs in 100 m bins
% now sum the traces in distance bins
distBin = min(rayStore(:,5)):20:max(rayStore(:,5));
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

figure(1);
hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
for i = 1:1:length(rayStore(:,1))
    rayX = [rayStore(i,1);rayStore(i,3)];
    rayY = [rayStore(i,2);rayStore(i,4)];
    
    plot(rayX,rayY,'k','LineWidth',2);
end
hold off;

%%
% get the fk transform
if(ifFK)
    newDistBin = (distBin(1,1:nBins-1) + distBin(1,2:nBins))/2;
    stackedTrace = stackedTrace(:,~isnan(stackedTrace(1,:)));
    interpDist = newDistBin;
    newDistBin = newDistBin(1,~isnan(stackedTrace(1,:)));
    
    interpStack = [];
    for i = 1:1:length(tArray)
        interpStack(i,:) = interp1(newDistBin,stackedTrace(i,:),interpDist);
    end
    % check if anty trailing NaN prssent
    interpStack = interpStack(:,~isnan(interpStack(1,:)));
    interpDist = interpDist(1,~isnan(interpStack(1,:)));
    
    [spec,f,kx]=fktran(interpStack,tArray,interpDist);
    % normalize every frequency bin by maximum
    fStart = 2; fEnd = 5;
    fStartInd = find(f>=fStart,1,'first');
    fEndInd = find(f>=fEnd,1,'first');
    fSmall = f(fStartInd:fEndInd,1);
    % for i = fStartInd:1:fEndInd
    %     spec(i,:) = abs(spec(i,:))/max(abs(spec(i,:)));
    % end
    spec = spec(fStartInd:fEndInd,:);
    figure(32)
    surf(kx,f(fStartInd:fEndInd),abs(spec));
    shading interp;
    view(2);
end