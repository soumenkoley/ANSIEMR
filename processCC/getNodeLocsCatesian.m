% this script was written to get the node locations in cartesian
% coordinates
clear; close all;

[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
stnLocs = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
refLong = 5.9060966; refLat = 50.7513927;

R = 6371000; % radius of earth in meters
for i = 1:1:length(allStn)
    stnInd = find(allStn==allStn(i,1));
    stnLocUse(i,1:3) = stnLocs(stnInd,1:3);
    stnLocUse(i,4) = stnInd;
    [xLong,yLat] = calculatedist(stnLocUse(i,1:2),[refLat,refLong],R);
    stnLocUse(i,5:6) = [xLong,yLat];
    %pairAz(i,3) = azimuth(pairNodeALatLong(1,1),pairNodeALatLong(1,2),pairNodeBLatLong(1,1),pairNodeBLatLong(1,2));
    
end

for i = 1:1:length(allStn)
    [xLong,yLat] = calculatedist(stnLocs(i,1:2),[refLat,refLong],R);
    nodeLocationsCartesian(i,1) = i;
    nodeLocationsCartesian(i,2:3) = [xLong,yLat]; 
    nodeLocationsCartesian(i,4) = stnLocs(i,3);
end

figure(1);
maxAlti = max(nodeLocationsCartesian(:,4));
minAlti = min(nodeLocationsCartesian(:,4));

maxAlti = maxAlti + 0.1*maxAlti;
minAlti = minAlti-0.1*minAlti;

colLev = 100;
colVal = jet(colLev);

figure(1);
hold on;
for i = 1:1:length(allStn)
    colInd = round((nodeLocationsCartesian(i,4)-minAlti)/(maxAlti-minAlti)*colLev);
    colNow = colVal(colInd,:);
    plot(nodeLocationsCartesian(i,2),nodeLocationsCartesian(i,3),'o','MarkerFaceColor',...
        colNow,'MarkerSize',10);
end
colorbar; colormap('jet');
caxis([minAlti,maxAlti]);
hold off;