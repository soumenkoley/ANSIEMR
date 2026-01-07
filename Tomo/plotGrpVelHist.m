% this script was written to plot the group velocity histograms as a
% function of frequency
clear; %close all;

fPathStore = 'B:\LimburgBigSurvey1CC-Pair\FinalGrpVelPics\';
load('grpTimesAll.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% get the not-NaN indices from the travelTime matrix
fSee = 2.4;

for fNo = 1:1:length(fSee)
    fSeeInd = find(fExt>=fSee(fNo),1,'first');

    notNaN = (~isnan(grpTimeStore(fSeeInd,:)));

    %%
    rayAttTemp = rayAttribute(notNaN,:);    
    grpTimeStoreTemp = grpTimeStore(fSeeInd,notNaN);
    grpTimeErrStoreTemp = grpTimeErrStore(fSeeInd,notNaN);
    grpSNRStoreTemp = grpSNRStore(fSeeInd,notNaN);

    velDist = rayAttTemp(:,5)./grpTimeStoreTemp';

    % make the plots and get the stats
    velMin = 500; velMax = 3000;
    velEdges = velMin:100:velMax;
    [N(:,fNo),x] = histcounts(velDist,velEdges);
    N(:,fNo) = N(:,fNo)/max(N(:,fNo));
end
velCent = (velEdges(1:(end-1)) + velEdges(2:end))/2;
% figure(1)
% subplot(1,2,1)
% surf(fSee,velCent,N);
% shading interp; view(2);
%%
% now do the analysis after SNR cut-off
for fNo = 1:1:length(fSee)
    fSeeInd = find(fExt>=fSee(fNo),1,'first');

    notNaN = (~isnan(grpTimeStore(fSeeInd,:)));

    %%
    rayAttTemp = rayAttribute(notNaN,:);    
    grpTimeStoreTemp = grpTimeStore(fSeeInd,notNaN);
    grpTimeErrStoreTemp = grpTimeErrStore(fSeeInd,notNaN);
    grpSNRStoreTemp = grpSNRStore(fSeeInd,notNaN);
    goodSNRInd = (grpSNRStoreTemp>=5);
    rayAttTemp = rayAttTemp(goodSNRInd,:);
    grpTimeStoreTemp = grpTimeStoreTemp(goodSNRInd);
    grpTimeErrStoreTemp = grpTimeErrStoreTemp(goodSNRInd);

    velDistNew = rayAttTemp(:,5)./grpTimeStoreTemp';

    % make the plots and get the stats
    [NNew(:,fNo),x] = histcounts(velDistNew,velEdges);
    NNew(:,fNo) = NNew(:,fNo)/max(NNew(:,fNo));
end
% figure(1);
% subplot(1,2,2)
% surf(fSee,velCent,NNew);
% shading interp; view(2);

