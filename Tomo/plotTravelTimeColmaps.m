% this script was written to observe the peak delay times with respect to
% some virtual source

clear;
close all;
% first load the cross-correlation file
load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\GFStore\AvgCorrAllDays.mat');
%load('D:\GreensFunction\psdCorrAvgDay4-15.mat');
%load('D:\GreensFunction\rayAttributeDay4-15.mat')
%load('D:\GreensFunction\lagValDay4-15.mat')
addpath('C:\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive2');
fSamp = 25;
dt = 1/fSamp;
fLowC = 2.5; fHighC=9;

refSensor = 87;
refSensorIndA = find(rayAttribute(:,1)==refSensor);
refSensorIndB = find(rayAttribute(:,2)==refSensor);
refSensorInd = [refSensorIndA;refSensorIndB];

good = 1;
for i = 1:1:length(refSensorInd)
    nodeCond = (rayAttribute(refSensorInd(i,1),2)==194) || ...
        (rayAttribute(refSensorInd(i,1),2)==94) ||(rayAttribute(refSensorInd(i,1),1)==194) ...
        || (rayAttribute(refSensorInd(i,1),1)==94) || (rayAttribute(refSensorInd(i,1),1)==1112)...
        || (rayAttribute(refSensorInd(i,1),2)==1112); 
    if(nodeCond)
        % do not store the data
        disp('I am here');
    else
        refSensorIndN(good,1) = refSensorInd(i,1);
        good = good+1;
    end
end
refSensorInd = refSensorIndN;

newRayAttribute = rayAttribute(refSensorInd,:);
ccData = psdCorrAvg(:,refSensorInd);
radialOffset = rayAttribute(refSensorInd,3);
tAxis = lagVal/fSamp;

% now only thing to decide is if these offsets are negative
% so just compare the y coordinate of points

% SNR parameters
vMin = 50; vMax = 400;
tN1 = 6; tN2 = 10;
snrThresh = 10;

cGood = 1; acGood = 1;
for i = 1:1:length(refSensorInd)
    n1 = newRayAttribute(i,1); n2 = newRayAttribute(i,2);
    distOut = newRayAttribute(i,3);
    if(n1==refSensor)
        % then the order is fine no need to do anything to ccData
        nodeLocPlot(i,:) = newRayAttribute(i,6:7);
    else
        ccData(:,i) = flipud(ccData(:,i));
        nodeLocPlot(i,:) = newRayAttribute(i,4:5);
    end
    [ccData(:,i)]=bandpass_n(ccData(:,i),fLowC,fHighC,1/fSamp,10);
    acausalData = [ccData(tAxis==0,i);flipud(ccData(tAxis<0,i))];
    causalData = ccData(tAxis>=0,i);
    
    SNRAcausal = SNR_cc(acausalData,vMin,vMax,distOut,tN1,tN2,dt);
    SNRCausal = SNR_cc(causalData,vMin,vMax,distOut,tN1,tN2,dt);
    % get the peak travel time, acausal first
    
    if(SNRAcausal>snrThresh)
        [acausalMag(acGood,1),acausalPeak(acGood,1)] = max(abs(acausalData));
        acausalPeak(acGood,1) = acausalPeak(acGood,1)/fSamp;
        nodeLocAcausal(acGood,:) = nodeLocPlot(i,:);
        acGood = acGood+1;
    end
    if(SNRCausal>snrThresh)
        [causalMag(cGood,1),causalPeak(cGood,1)] = max(abs(causalData));
        causalPeak(cGood,1)=causalPeak(cGood,1)/fSamp;
        nodeLocCausal(cGood,:) = nodeLocPlot(i,:);
        cGood = cGood+1;
    end
        
end

% now define colorbar
colLevels = 500;
colMat = jet(colLevels);
%maxT = max(acausalPeak)+0.2; minT = min(acausalPeak)-0.1;
maxT = 4; minT = 0.02;
figure(1)
subplot(1,2,1)
% plot for acausal times
hold on;
for i = 1:1:length(acausalPeak)
    colInd = floor(((acausalPeak(i,1)-minT)/(maxT-minT))*colLevels);
    colInd = min(colInd,colLevels);
    colNow = colMat(colInd,:);
    plot(nodeLocAcausal(i,1),nodeLocAcausal(i,2),'o','color',colNow,'MarkerFaceColor',colNow,'MarkerSize',12);
end
colorbar;
colormap('jet'); caxis([minT,maxT]);
hold off;

subplot(1,2,2)
% plot for acausal times
%maxT = max(causalPeak)+0.1; minT = min(causalPeak)-0.1;
maxT = 4; minT = 0.02;

hold on;
for i = 1:1:length(causalPeak)
    colInd = floor(((causalPeak(i,1)-minT)/(maxT-minT))*colLevels);
    colInd = min(colInd,colLevels);
    colNow = colMat(colInd,:);
    plot(nodeLocCausal(i,1),nodeLocCausal(i,2),'o','color',colNow,'MarkerFaceColor',colNow,'MarkerSize',12);
end
colorbar;
colormap('jet'); caxis([minT,maxT]);
hold off;

% plotting corelation magnitudes
% now define colorbar
colLevels = 500;
colMat = jet(colLevels);
maxAmp = max(acausalMag)+0.005; minAmp = min(acausalMag)-0.005;
%maxT = 4; minT = 0.02;
figure(2)
subplot(1,2,1)
% plot for acausal times
hold on;
for i = 1:1:length(acausalPeak)
    colInd = floor(((acausalMag(i,1)-minAmp)/(maxAmp-minAmp))*colLevels);
    colNow = colMat(colInd,:);
    plot(nodeLocAcausal(i,1),nodeLocAcausal(i,2),'o','color',colNow,'MarkerFaceColor',colNow,'MarkerSize',12);
end
colorbar;
colormap('jet'); caxis([minAmp,maxAmp]);
hold off;

subplot(1,2,2)
% plot for acausal times
maxAmp = max(causalMag)+0.005; minT = min(causalMag)-0.005;
hold on;
for i = 1:1:length(causalPeak)
    colInd = floor(((causalMag(i,1)-minAmp)/(maxAmp-minAmp))*colLevels);
    colNow = colMat(colInd,:);
    plot(nodeLocCausal(i,1),nodeLocCausal(i,2),'o','color',colNow,'MarkerFaceColor',colNow,'MarkerSize',12);
end
colorbar;
colormap('jet'); caxis([minAmp,maxAmp]);
hold off;
