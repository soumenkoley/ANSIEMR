% this B:\LimburgBigSurvey1CC-Pair\SubArrayCCscript was written to get the theoretical phase travel times between
% station pairs based on a group or phase velocity map at a particular
% frequency.
% the phase/group velocity maps were derived previously

clear;
close all;
%fPath = 'B:\LimburgBigSurvey1CC-Pair\velMaps\';
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelStarsApril.mat';
fExtN = 1.0:0.05:5.0;
% for i = 1:1:length(fExt)
%     fName = ['grpVel',num2str(fExt(i)),'.mat'];
%     load([fPath,fName]);
%     mSlow(:,i) = mFinalNew;
% end
% load the theoretical group velocity maps
load(fPath);
% grpMat has dimensions nBox x nFreq (1122 x 75 in this case)
mSlow = 1./grpMat;
% however fExt and fExtN might not be the same, so align them
fStartInd = find(fExt>=fExtN(1,1),1,'first');
fEndInd = find(fExt>=fExtN(1,end),1,'first');
mSlow = mSlow(:,fStartInd:fEndInd);
% load the allCCs file
load(['B:\LimburgBigSurvey1CC-Pair\SubArrayCC\','allCCFundUpdApril.mat']);
dx = 200; dy = 200;
minRayLength = 1; % units in meters

% set up the new coordinate axis
[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);

figure(1)
subplot(1,2,1)
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;

[s1,s2] = size(rayAttribute);
nBoxes = (floor(xLimit/dx))*(floor(yLimit/dy));
rayCount = zeros(nBoxes,1);
for rayNo = 1:1:s1
    nodeA = rayAttribute(rayNo,7); nodeB = rayAttribute(rayNo,8);
    nodeAInd = find(newLoc(:,1)==nodeA); nodeBInd = find(newLoc(:,1)==nodeB);
    nodeALoc = newLoc(nodeAInd,2:3); nodeBLoc = newLoc(nodeBInd,2:3);
    rayStart = [nodeALoc(1,1),nodeBLoc(1,1)]';
    rayEnd = [nodeALoc(1,2),nodeBLoc(1,2)]';
    %plot(rayStart,rayEnd,'k','LineWidth',1);
    [propDetails] = findTravelBoxes(nodeALoc,nodeBLoc,dx,dy,xLimit,yLimit);
    rayCount(propDetails(:,4),1) = rayCount(propDetails(:,4),1)+1;
    propDetails(propDetails(:,5)<minRayLength,5) = 0;
    kernelMat = zeros(1,nBoxes);
    kernelMat(1,propDetails(:,4)) = propDetails(:,5);
    % now find the theoretical travel time for each ray
    for f = 1:1:length(fExtN)
        thTime(f,rayNo) = kernelMat*mSlow(:,f);
    end
    %plot(propDetails(:,1),propDetails(:,2),'ro','MarkerSize',4,'MarkerFaceColor','r');
    %disp('Doing one ray!');
end
hold off;

for i = 1:1:s1
    thVelNow(:,i) = rayAttribute(i,5)./thTime(:,i);
end
figure(2);
plot(fExt,thVelNow);
