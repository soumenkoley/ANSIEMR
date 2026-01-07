% this script was written to assemble the group velocity picks and prepare
% it for Tomography
clear; close all;

fPathStore = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\HighVel\AllFreqPicksPerc30\';
%load('grpTimesAllNew.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\HighVel\GrpTimePickAllPerc30.mat')
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% get the not-NaN indices from the travelTime matrix
fSee = 2;
fStoreName = ['grpTimeNew',num2str(fSee),'.mat'];
fSeeInd = find(fExt>=fSee,1,'first');

notNaN = (~isnan(grpTimeStore(fSeeInd,:)));

%%
rayAttTemp = rayAttribute(notNaN,:);
grpTimeStoreTemp = grpTimeStore(fSeeInd,notNaN);
grpTimeErrStoreTemp = grpTimeErrStore(fSeeInd,notNaN);
grpSNRStoreTemp = grpSNRStore(fSeeInd,notNaN);

velDist = rayAttTemp(:,5)./grpTimeStoreTemp';

% make the plots and get the stats
velMin = min(velDist)-100; velMax = max(velDist)+100;

nCol = 100;
colVal = jet(nCol);
figure(1)
subplot(1,2,1)
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ko','MarkerFaceColor','k');

hold on;
for i = 1:1:length(velDist)
    velNow = velDist(i,1);
    velInd = floor((velNow-velMin)/(velMax-velMin)*nCol);
    colNow = colVal(velInd,:);
    rayX = [rayAttTemp(i,1),rayAttTemp(i,3)];
    rayY = [rayAttTemp(i,2),rayAttTemp(i,4)];
    plot(rayX,rayY,'color',colNow);
end
hold off;
colorbar;colormap('jet');
caxis([velMin,velMax]);

figure(2)
subplot(1,2,1)
velEdges = velMin:100:velMax;
histogram(velDist,velEdges);
%%
% now do the analysis after SNR cut-off
goodSNRInd = (grpSNRStoreTemp>=5);
rayAttTemp = rayAttTemp(goodSNRInd,:);
grpTimeStoreTemp = grpTimeStoreTemp(goodSNRInd);
grpTimeErrStoreTemp = grpTimeErrStoreTemp(goodSNRInd);

velDistNew = rayAttTemp(:,5)./grpTimeStoreTemp';

% now do the plots to check how everything works

velMin = min(velDistNew)-100; velMax = max(velDistNew)+100;

figure(1)
subplot(1,2,2)
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ko','MarkerFaceColor','k');

hold on;
for i = 1:1:length(velDistNew)
    velNow = velDistNew(i,1);
    velInd = floor((velNow-velMin)/(velMax-velMin)*nCol);
    colNow = colVal(velInd,:);
    rayX = [rayAttTemp(i,1),rayAttTemp(i,3)];
    rayY = [rayAttTemp(i,2),rayAttTemp(i,4)];
    plot(rayX,rayY,'color',colNow);
end
hold off;
colorbar;colormap('jet');
caxis([velMin,velMax]);

figure(2)
subplot(1,2,2)
velEdges = velMin:100:velMax;
histogram(velDistNew,velEdges);

%% prepare to store for usage in Tomography
for i = 1:1:length(rayAttTemp(:,1))
    rayAttributeStoreFull(1,i) = rayAttTemp(i,7);
    rayAttributeStoreFull(2,i) = rayAttTemp(i,8);
    rayAttributeStoreFull(3,i) = rayAttTemp(i,5);
    rayAttributeStoreFull(4,i) = rayAttTemp(i,1);
    rayAttributeStoreFull(5,i) = rayAttTemp(i,2);
    rayAttributeStoreFull(6,i) = rayAttTemp(i,3);
    rayAttributeStoreFull(7,i) = rayAttTemp(i,4);
    rayAttributeStoreFull(8,i) = rayAttTemp(i,6);
end
phaseVelStoreFull = velDistNew';
timePickStoreFull = grpTimeStoreTemp;
errorTStoreFull = grpTimeErrStoreTemp;

fStoreFull = [fPathStore,fStoreName];
save(fStoreFull,'rayAttributeStoreFull','phaseVelStoreFull',...
    'timePickStoreFull','errorTStoreFull');