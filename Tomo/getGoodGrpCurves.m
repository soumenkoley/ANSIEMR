% this script was written to check which station pairs are good
% implying they have at leasy 40 good points with SNR >=7
% we will also check which one of them have at least 20 continuous picks
% station separation at least 1 km
clear; close all;

%% load input files
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');

% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');

% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% load the GV picks
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\LowVelAug\GrpTimePick.mat');

fAll = 0.5:0.05:5;
fSee = 1:0.05:5;
fStartInd = find(fAll>=fSee(1),1,'first');
fEndInd = find(fAll>=fSee(end),1,'first');
grpSNRStore = grpSNRStore(fStartInd:fEndInd,:);
nRays = size(grpSNRStore,2);
% now loop and find how many rays have at least minDispPoints
minDist = 1000; % one km
minDispPoints = 40;
snrCut = 7;

grpSNRStore = grpSNRStore(:,rayAttribute(:,5)>=minDist);
grpTimeStore = grpTimeStore(:,rayAttribute(:,5)>=minDist);
grpTimeErrStore = grpTimeErrStore(:,rayAttribute(:,5)>=minDist);
rayAttribute = rayAttribute(rayAttribute(:,5)>=minDist,:);

snr7Count = sum(grpSNRStore>snrCut,1,'omitnan');
rayAttGood = rayAttribute(snr7Count>=minDispPoints,:);

% now actually get those station names and pairs
figure(1);
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),...
    'bo','MarkerSize',6,'MarkerFaceColor','b');
hold on;
for i = 1:1:length(rayAttGood)
    rayX = [rayAttGood(i,1);rayAttGood(i,3)];
    rayY = [rayAttGood(i,2);rayAttGood(i,4)];
    plot(rayX,rayY,'k');
    stnPairs(i,1:2) = [allStn(rayAttGood(i,7),1),allStn(rayAttGood(i,8),1)];
end
hold off;

