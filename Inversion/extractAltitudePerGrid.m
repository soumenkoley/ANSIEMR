% this script was written to get the altitude per grid point
clear; close all;

% load the node locations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% convert to new coordinate system since the tomography is done on the new
% grid map
minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newLoc(:,1) = nodeLocationsCartesian(:,1);
newLoc(:,2) = nodeLocationsCartesian(:,2)-minX;
newLoc(:,3) = nodeLocationsCartesian(:,3)-minY;

% load the altitude val
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\altitudeVal.mat');

% load the XYGrid points
load('AllXYGridPoints.mat');

nX = length(xCoord);
nY = length(yCoord);
% the model matrix will be (nY x nX) matrix
% travel across x, then y to create all set of points

gC = 1;
for i = 1:1:nY
    for j = 1:1:nX
        xyAllGrid(gC,1) = xCoord(j);
        xyAllGrid(gC,2) = yCoord(i);
        gC = gC+1;
    end
end

% now interpolate 2D
altIntp = griddata(newLoc(:,2),newLoc(:,3),altiVal,xyAllGrid(:,1),xyAllGrid(:,2));


