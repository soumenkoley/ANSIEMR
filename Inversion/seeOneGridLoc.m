% this sript was just written as a test to see the location of a grid point
% based on the grid number, and then also see the location of the stations
% basically you can select the correct grid point based on observation
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\InfoStore\nodeLocationsCartesian.mat');

minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newLoc(:,1) = nodeLocationsCartesian(:,1);
newLoc(:,2) = nodeLocationsCartesian(:,2)-minX;
newLoc(:,3) = nodeLocationsCartesian(:,3)-minY;

maxX = max(newLoc(:,2));
maxY = max(newLoc(:,3));

figure(1)
hold on;
for i = 1:1:length(newLoc)
    plot(newLoc(i,2),newLoc(i,3),'bo');
    %text(newLoc(i,2),newLoc(i,3),['--',num2str(newLoc(i,1))]);
end
%hold off;

% now plot the grid location
gridNo = 2680;

gridPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey2\assembleCrossCorr\Tomo\doInver\SNR5\';
load([gridPath,'grid',num2str(gridNo),'\','grpVel.mat']);

plot(xyNow(1,1),xyNow(1,2),'ko','MarkerSize',6,'MarkerFaceColor','k');
