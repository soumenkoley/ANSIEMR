% this script was written tp generate the theoretical group velocity
% maps;We use the group velocity pcked per star, generate the grid, and
% use Gaussian smoothing to generate the velocity maps at each frequency

clear; close all;

% load the node locations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

fExt = 1.0:0.05:5.0;
nFreq = length(fExt');
fVec = 1.0:0.05:5.0;
fVec = fVec';
% now work with the stars
%starVal = [1:1:23,25,26,29:34];
%starVal = 1:1:14;
%starVal = 1:1:46; %redoing it for the first StarsJan case, picks were good
%starHigh = [1:6,8:10,35,27:29];
%starVal = setdiff(starVal,starHigh);

% use the below star val for FundFeb
%starVal = [7,11:26,30:34,36:43,47,48];

% use the below star vals for FundUpdApril
starVal = [1:43,45:51];


% use the below starval for HighFeb
%starVal = [1:7,9:13];
nStars = length(starVal);

% convert the star locations from lat long to cartesian, sing the same
% reference lat long that was used to generate nodeLocationsCartesian.mat
refLong = 5.9060966; refLat = 50.7513927;

R = 6371000; % radius of earth in meters
%starPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\StarsJan\FundFeb\';
starPath = ['A:\TestInver\StarsApril\'];

for i = 1:1:nStars
    fPath = [starPath,'Star',num2str(starVal(i)),'\starLocs.txt'];
    starLoc = load(fPath);
    
    [xLong,yLat] = calculatedist([starLoc(1,1),starLoc(2,1)],[refLat,refLong],R);
    starLocsCartesian(i,1) = i;
    starLocsCartesian(i,2:3) = [xLong,yLat];
    %pairAz(i,3) = azimuth(pairNodeALatLong(1,1),pairNodeALatLong(1,2),pairNodeBLatLong(1,1),pairNodeBLatLong(1,2));
end

figure(1); subplot(1,2,1);hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo','MarkerSize',6);
plot(starLocsCartesian(:,2),starLocsCartesian(:,3),'ro','MarkerSize',6,'MarkerFaceColor','r');
hold off;

%% set up the grid
dx = 200; dy = 200;
sigmaX = 200; sigmaY = 200; % units in meters
[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);
nBoxes = (floor(xLimit/dx))*(floor(yLimit/dy));

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

xCoord = ((1:1:nBoxRow)-1)*dx + dx/2;
yCoord = ((1:1:nBoxCol)-1)*dy + dy/2;

% also time to transform the locations of the stars

minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newStarLocs(:,1) = starLocsCartesian(:,1);
newStarLocs(:,2) = starLocsCartesian(:,2)-minX;
newStarLocs(:,3) = starLocsCartesian(:,3)-minY;

% plot again after coordinate transformation
figure(1); subplot(1,2,2);hold on;
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerSize',6);
plot(newStarLocs(:,2),newStarLocs(:,3),'ro','MarkerSize',6,'MarkerFaceColor','r');
hold off;

%% now perform the Gaussian smoothing to prepare theoretical group velocity maps

wMat = zeros(nBoxes,length(fVec));
grpMat = zeros(nBoxes,length(fVec));

figure(3); hold on;
for i = 1:1:nStars
    fPath = [starPath,'Star',num2str(starVal(i)),'\GrpUpdAprilTheo.txt'];
    grpVel = load(fPath);
    % interpolate grpVel on fVec
    grpVelN = interp1(grpVel(:,1),grpVel(:,2),fVec);
    grpVel = [fVec,grpVelN];
    fStartInd = find(grpVel(:,1)>=fExt(1,1),1,'first');
    fEndInd = find(grpVel(:,1)>=fExt(end),1,'first');
    grpVelUse = grpVel(fStartInd:fEndInd,1:2);
    expVal = (boxCoord(:,1)-newStarLocs(i,2)).^2/(2*sigmaX^2) + ...
        (boxCoord(:,2)-newStarLocs(i,3)).^2/(2*sigmaY^2);
    plot(grpVelUse(:,1),grpVelUse(:,2),'b','LineWidth',2);
    gauss2D = exp(-expVal);
%     figure(3);
%     surf(xCoord,yCoord,vec2mat(gauss2D,nBoxRow));
%     view(2);
    wMat = wMat + repmat(gauss2D,1,length(fVec));
    grpMat = grpMat + repmat(grpVelUse(:,2)',nBoxes,1).*repmat(gauss2D,1,length(fVec));
    
    %disp('Stop!');
end
hold off;

grpMat = grpMat./wMat;

fSee = 2.0; % units in Hz
fSeeInd = find(fExt>=fSee,1,'first');
velMat = vec2mat(grpMat(:,fSeeInd),nBoxRow);
figure(289)
hold on;
imagesc(xCoord,yCoord,velMat);
%shading interp; view(2);
colorbar; colormap('jet');
%caxis([100,3500]);
plot(newLoc(:,2),newLoc(:,3),'ko','MarkerFaceColor','k','MarkerSize',6);
hold off;