% this script was written to see the virtual noise gather
clear; close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive2');

load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\ccAssembled\RefStationWUJ2A.mat');
load('HdBp0_2to0_8Hz.mat');
A = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
refStn = "WUJ2A";
refStnInd = find(allStn==refStn,1,'first');
refLatLong = A(refStnInd,1:2);

R = 6371000; % in meters
for i = 1:1:length(endStns)
    endStnInd = find(allStn==endStns(i,1),1,'first');
    latLongNow = A(endStnInd,1:2);
    [xlong,ylat] = calculatedist(latLongNow,refLatLong,R);
    nodeLocCart(i,1:2) = [xlong,ylat];
    nodeLocCart(i,3) = sqrt(xlong^2 + ylat^2);
end

load('HdBp0_2to0_8Hz.mat');
%ccAvg = filtfilt(HdBp0_2to0_8Hz.Numerator,1,ccAvg);

for i = 1:1:length(ccAvg(1,:))
    ccAvg(:,i) = filtfilt(HdBp0_2to0_8Hz.Numerator,1,ccAvg(:,i));
    ccAvg(:,1) = ccAvg(:,i)/max(abs(ccAvg(100:900,i)));
end

fSamp = 25;
tArray = [(-400:1:-1),(0:1:400)]/fSamp;
plotseis(ccAvg(100:900,:),tArray,nodeLocCart(:,3),1);

% this lat long will be used to plot the sensors on the map
refLong = 5.9060966; refLat = 50.7513927;
tArrayNew = [(-300:1:-1),(0:1:300)]/fSamp;
%tArrayNew = [(-300:1:-1)]/fSamp;
%tArrayNew = [0:1:399]/fSamp;
for i = 1:1:length(endStns)
    endStnInd = find(allStn==endStns(i,1),1,'first');
    latLongNow = A(endStnInd,1:2);
    [xlongNew,ylatNew] = calculatedist(latLongNow,[refLat,refLong],R);
    nodeLocCartNew(i,1:2) = [xlongNew,ylatNew];
    % now get the peak travel time
    %ccAvgNow = [ccAvg(501,i);(flipud(ccAvg(1:1:500,i)) + ccAvg(502:1:1001,i))];
    %[M,I] = max(abs(ccAvg(200:499,i)));
    [M,I] = max(abs(ccAvg(200:800,i)));
    if(ccAvg(1,i)~=0)
        tMax(i,1) = tArrayNew(I);
    end
end

minT = min(tMax)-0.1; maxT = max(tMax)+0.1;
figure(5)
hold on;
nColors = 100;
colLev = jet(nColors);
for i = 1:1:length(tMax)
    colInd(i,1) = floor((tMax(i,1)-minT)/(maxT-minT)*nColors);
    plot(nodeLocCartNew(i,1),nodeLocCartNew(i,2),'o','MarkerSize',14,...
        'MarkerFaceColor',colLev(colInd(i,1),:));
end
hold off;
colormap('jet'); colorbar;
caxis([minT,maxT]);


