function [newLoc,xLimit,yLimit] = setUpGeometry(dx,dy)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% figure(1)
% hold on;
% for i = 1:1:length(goodNodeLocations)
%     plot(goodNodeLocations(i,2),goodNodeLocations(i,3),'bo');
%     text(goodNodeLocations(i,2),goodNodeLocations(i,3),['--',num2str(goodNodeLocations(i,1))]);
% end
% hold off;

minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newLoc(:,1) = nodeLocationsCartesian(:,1);
newLoc(:,2) = nodeLocationsCartesian(:,2)-minX;
newLoc(:,3) = nodeLocationsCartesian(:,3)-minY;

maxX = max(newLoc(:,2));
maxY = max(newLoc(:,3));

xLimit = (ceil(maxX/dx))*dx;
yLimit = (ceil(maxY/dy))*dy;
% figure(2)
% hold on;
% for i = 1:1:length(goodNodeLocations)
%     plot(newLoc(i,2),newLoc(i,3),'bo');
%     text(newLoc(i,2),newLoc(i,3),['--',num2str(newLoc(i,1))]);
% end
% hold off;

end

