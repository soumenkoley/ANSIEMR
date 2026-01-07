% this script computes the rayCount at any particular frequency
clear;
close all;
% load the rayStore file
load('D:\GreensFunction\GroupVelocityPicks\GroupVelocityPickDec27.mat');

dx = 25; dy = 25;
fGet = 2.9; % units in Hz
fExt = 2:0.1:8;
fGetInd = find(fExt >= fGet,1,'first');
[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);

figure(1)
subplot(1,2,1)
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;

[s1,s2] = size(rayAttributeStore);
nBoxes = (floor(xLimit/dx))*(floor(yLimit/dy));
rayCount = zeros(nBoxes,1);
for rayNo = 1:1:s1
    if(~isnan(groupVelStore(fGetInd,rayNo)))
        nodeA = rayAttributeStore(rayNo,1); nodeB = rayAttributeStore(rayNo,2);
        nodeAInd = find(newLoc(:,1)==nodeA); nodeBInd = find(newLoc(:,1)==nodeB);
        nodeALoc = newLoc(nodeAInd,2:3); nodeBLoc = newLoc(nodeBInd,2:3);
        rayStart = [nodeALoc(1,1),nodeBLoc(1,1)]';
        rayEnd = [nodeALoc(1,2),nodeBLoc(1,2)]';
        plot(rayStart,rayEnd,'k','LineWidth',2);
        [propDetails] = findTravelBoxes(nodeALoc,nodeBLoc,dx,dy,xLimit,yLimit);
        rayCount(propDetails(:,4),1) = rayCount(propDetails(:,4),1)+1;
        %plot(propDetails(:,1),propDetails(:,2),'ro','MarkerSize',4,'MarkerFaceColor','r');
    end
end
hold off;

% now get the coordinates of the boxes
boxCoord = zeros(nBoxes,2);
for boxNo = 1:1:nBoxes
    nBoxRow = floor(xLimit/dx);
    nBoxCol = floor(yLimit/dy);
    if(rem(boxNo,nBoxRow)==0)
        colNo = nBoxRow;
        rowNo = floor(boxNo/nBoxRow);
    else
        colNo = rem(boxNo,nBoxRow);
        rowNo = floor(boxNo/nBoxRow)+1;
    end
    boxCoord(boxNo,1) = (colNo-1)*dx + dx/2;
    boxCoord(boxNo,2) = (rowNo-1)*dy + dy/2;
end

rayCountMat = vec2mat(rayCount',nBoxRow);
xCoord = ((1:1:nBoxRow)-1)*dx + dx/2;
yCoord = ((1:1:nBoxCol)-1)*dy + dy/2;

figure(1)
subplot(1,2,2)
surf(xCoord,yCoord,rayCountMat);
caxis([0,10]);
colormap('jet');
shading interp;
