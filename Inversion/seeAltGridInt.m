% this script was written to see the altitude variation for a particular
% grid interval
clear; close all;
gridStart = 855; gridEnd = 877;

% load the altitude values
load('AllGridAlti');
% load the All XY grid points
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

gC = 1;
for i = gridStart:gridEnd
    altNow(gC,1) = altIntp(i,1);
    xDist(gC,1) = xyAllGrid(i,1);
    yDist(gC,1) = xyAllGrid(i,2);
    gC = gC+1;
end

figure(1);
subplot(1,2,1);
load('newNodeLocs.mat');
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;
plot(xDist,yDist,'k*');
hold off;

subplot(1,2,2)
plot(xDist,altNow,'b','LineWidth',2);
ylim([0,400])
