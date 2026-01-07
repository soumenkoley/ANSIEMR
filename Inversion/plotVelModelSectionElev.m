% thi script was written to display he velocity model as a 2D
% along a line passing through some grid points
% elevation is also taken care of

clear; close all;
% load the altitude information
load('AllGridAlti.mat');

gridStart = 923; gridEnd = 944;
NBest = 500;
% extract the alti information
gC = 1;
for i = gridStart:1:gridEnd
    altNow(gC,1) = altIntp(i,1);
    if(isnan(altNow(gC,1)))
        altNow(gC,1) = altNow(gC-1,1);
    end
    gC = gC+1;
end

dx = 200; dz = 1; % units in meters
dzMax = 900;
dzArray = 0:dz:dzMax;
maxAlti = 400;
gridVal = gridStart:1:gridEnd;
nGrid = length(gridVal);
dxArray = 0:dx:((nGrid-1)*dx);

gPath = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvertApril\grid';

velMat = zeros(length(dzArray),nGrid);
velMatMean = zeros(length(dzArray),nGrid);
for i = 1:1:nGrid
    % load the inversion output
    fPathFull = [gPath,num2str(gridVal(i)),'\invOutApril.mat'];
    if(exist(fPathFull))
        load([gPath,num2str(gridVal(i)),'\invOutApril.mat']);
        
        % time to get the mean model
        %vsMean = mean(vs,2);
        for j = 1:1:length(dzArray)
            zInd = find(h(:,1)>=dzArray(j),1,'first');
            velMat(j,i) = vs(zInd,1);
            % now loop across all the 500 models
            for modNo = 1:1:NBest
                zInd = find(h(:,modNo)>=dzArray(j),1,'first');
                velMatMean(j,i) = velMatMean(j,i) + vs(zInd,modNo);
            end
            
        end
        velMatMean(:,i) = velMatMean(:,i)/NBest;
        %     p = polyfit(dzArray',velMat(:,i),3);
        %     velMatPolyFit(:,i) = p(1)*dzArray.^3 + p(2)*dzArray.^2+ p(3)*dzArray + p(4);
        %     disp('one grid');
        relErrAll(i,1) = min(relErr);
    else
        velMat(:,i) = velMat(:,i-1);
        velMatMean(:,i) = velMatMean(:,i-1);
        relErrAll(i,1) = relErrAll(i-1,1);
    end
end
% for i = 1:1:length(dzArray)
%     velMat(i,:) = smooth(velMat(i,:),5);
%     velMatMean(i,:) = smooth(velMatMean(i,:),5);
% end

% make the dzArray negative
maxdZ = max(dzArray);
dzArray = -dzArray;
newdZArray = maxAlti:-1:-500;
velMatElev = NaN*ones(length(newdZArray),nGrid);
velMatMeanElev = NaN*ones(length(newdZArray),nGrid);

for i = 1:1:length(gridVal)
    dzArrayNow = dzArray+altNow(i);
    maxdZInd = find(dzArrayNow<=-500,1,'first');
    for j = 1:1:maxdZInd
        ddInd = find(newdZArray<=dzArrayNow(1,j),1,'first');
        velMatElev(ddInd,i) = velMat(j,i);
        velMatMeanElev(ddInd,i) = velMatMean(j,i);
    end
    %disp('Stop');
end

for i = 1:1:length(dzArray)
    velMat(i,:) = mySmooth(velMat(i,:),3);
    velMatMean(i,:) = mySmooth(velMatMean(i,:),3);
    velMatElev(i,:) = mySmooth(velMatElev(i,:),3);
    velMatMeanElev(i,:) = mySmooth(velMatMeanElev(i,:),3);
end

figure(1);
subplot(1,2,1)
surf(dxArray,dzArray,velMat);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');

figure(1);
subplot(1,2,2)
surf(dxArray,dzArray,velMatMean);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');

figure(2);
subplot(1,2,1)
surf(dxArray,newdZArray,velMatElev);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');

figure(2);
subplot(1,2,2)
surf(dxArray,newdZArray,velMatMeanElev);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');
figure(3);
plot(dxArray,relErrAll);

%%
% plot the grid locations
load('newNodeLocs.mat');
figure(4);
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerSize',6,'MarkerFaceColor','b');
hold on;
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

for i = gridStart:gridEnd
    plot(xyAllGrid(i,1),xyAllGrid(i,2),'k*');
end
hold off;
%% Observe velocities at depth
depth1 = -100; depth2 = -200; depth3 = -300;

depth1Ind = find(newdZArray<=depth1,1,'first');
depth2Ind = find(newdZArray<=depth2,1,'first');
depth3Ind = find(newdZArray<=depth3,1,'first');

velVecDepth1 = velMatElev(depth1Ind,:);
velVecDepth2 = velMatElev(depth2Ind,:);
velVecDepth3 = velMatElev(depth3Ind,:);

figure(5); hold on;
plot(dxArray,velVecDepth1,'b','LineWidth',2);
plot(dxArray,velVecDepth2,'r','LineWidth',2);
plot(dxArray,velVecDepth3,'k','LineWidth',2);
hold off;