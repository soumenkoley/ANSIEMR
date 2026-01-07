% this script was written to get a depth slice of the velocity model
clear; close all;
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\doInver\';
gridStartVal = [181,196,215,248,281,315,349,380,417,449,482,518,549,585,...
    622,657,692,726,759,793,827,835];
gridEndVal = [186,201,238,271,304,337,374,407,440,474,507,540,575,607,...
    641,672,705,740,774,804,830,838];

seeDepth = 100;
depthFromSurf = 300;

gridA = 585; gridB = 607;
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

velSliceVec = NaN*ones(nY*nX,1);
velSliceSurfVec = NaN*ones(nY*nX,1);

nGrid = length(gridStartVal);

for i = 1:1:nGrid
    fName = [fPath,'grid',num2str(gridStartVal(i)),'-',num2str(gridEndVal(i)),'\','velOut.mat'];
    load(fName);
    gridInterval = gridStartVal(i):gridEndVal(i);
    for j = 1:1:length(gridInterval)
        zInd = find(velStruct.newdZArray<=seeDepth,1,'first');
        velSliceVec(gridInterval(j),1) = velStruct.velMatElev(zInd,j);
        firstNotNanInd = find(~isnan(velStruct.velMatElev(:,j)), 1, 'first');
        velSliceSurfVec(gridInterval(j),1) = velStruct.velMatElev(firstNotNanInd +depthFromSurf,j);
    end
    disp('one done!');
    
end

% convert the vector to matrix

velSliceMat = vec2mat(velSliceVec,nX);
velSliceSurfMat = vec2mat(velSliceSurfVec,nX);

figure(1);
h = imagesc(xCoord,yCoord,velSliceMat);
set(h, 'AlphaData', ~isnan(velSliceMat))
shading interp;
colormap('jet');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
set(gca,'YDir','normal');
grid on; box on;

figure(2);
h = imagesc(xCoord,yCoord,velSliceSurfMat);
set(h, 'AlphaData', ~isnan(velSliceSurfMat))
shading interp;
colormap('jet');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
set(gca,'YDir','normal');
grid on; box on;

%%
% create the station location array
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\allStn.mat');
load('nodeLocTomo.mat');
 
 notPlotList = ["WKQQA";"W3QZA";"W2RAA";"XFRFA";"W4QGA";"XPQJA";...
    "W7PSA";"W4L7A";"XAL4A";"XELYA";"W9MRA";"XFMHA";"XMMBA";...
    "XFMWA"];
% add plot or not plot thing in the nodeNames cellarray

subArrayList = ["W1KLA";"W1LKA";"W3JWA";"W8KFA";"WUJ2A";"X0K9A";"X2KQA";...
    "XBKUA";"XDLDA";"XHKAA";"XUKKA";"XUKWA";"XYKJA";"Y0KWA";"Y8KWA";...
    "YDK3A";"YGKGA";"YHKVA";"ZDK2A";"ZOKUA";"ZUJ5A"];
pCount = 1;
for i = 1:1:length(allStn)
    ss = find(notPlotList == allStn{i,1});
    if(isempty(ss))
        newLocPlot(pCount,1:2) = newLoc(i,2:3);
        pCount = pCount+1;
    end
end
figure(1); hold on;
plot3(newLocPlot(:,1),newLocPlot(:,2),10000*ones(length(newLocPlot),1),'ko',...
    'MarkerSize',4,'MarkerFaceColor','k');
hold off;

figure(2); hold on;
plot3(newLocPlot(:,1),newLocPlot(:,2),10000*ones(length(newLocPlot),1),'ko',...
    'MarkerSize',4,'MarkerFaceColor','k');
hold off;
%%
%plot the vertical vel section between gridA and gridB
fName = [fPath,'grid',num2str(gridA),'-',num2str(gridB),'\','velOut.mat'];
load(fName);
dxArray = xyAllGrid(gridA:gridB,1);
figure(3);
surf(dxArray,velStruct.newdZArray,velStruct.velMatElev);
shading interp;
view(2);
set(gca,'YDir','normal');
colorbar;colormap('hsv');

%% plot the vel section on the depth slice
figure(4);
h = imagesc(xCoord,yCoord,velSliceSurfMat);
set(h, 'AlphaData', ~isnan(velSliceSurfMat))
shading interp;
colormap('jet');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
set(gca,'YDir','normal');
grid on; box on;
hold on;
plot3(xyAllGrid(gridA:gridB,1),xyAllGrid(gridA:gridB,2),...
    10000*ones(gridB-gridA+1,1),'k','LineWidth',2);
hold off;