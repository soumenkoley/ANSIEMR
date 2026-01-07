% this sript was just written as a test to see the location of a grid point
% based on the grid number, and then also see the location of the stations
% basically you can select the correct grid point based on observation
clear; close all

load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newLoc(:,1) = nodeLocationsCartesian(:,1);
newLoc(:,2) = nodeLocationsCartesian(:,2)-minX;
newLoc(:,3) = nodeLocationsCartesian(:,3)-minY;

maxX = max(newLoc(:,2));
maxY = max(newLoc(:,3));

figure(178)
hold on;
for i = 1:1:length(newLoc)
    plot(newLoc(i,2),newLoc(i,3),'bo');
    %text(newLoc(i,2),newLoc(i,3),['--',num2str(newLoc(i,1))]);
end
hold off;

fig1 = openfig('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\GrpTimePickFigs\FundVelApril1_6Hz.fig');
% now plot the grid location
gridNo = 1059;

for i = 1:1:length(gridNo)
    gridPath = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvertApril\';
    load([gridPath,'grid',num2str(gridNo(i)),'\','grpVel.mat']);
    phVel = load([gridPath,'grid',num2str(gridNo(i)),'\','phVelNowSmooth.txt']);
    figure(fig1); hold on;
    plot3(xyNow(1,1),xyNow(1,2),10000,'ko','MarkerSize',8);
    hold off;
    
    figure(179); hold on;
    plot(1.0:0.05:5,velNow(:,1),'b','LineWidth',2);
    hold off;
    
    figure(180); hold on;
    plot(phVel(:,1),phVel(:,2),'b','LineWidth',2);
    hold off;
    
end
