% this script was written to assemble the group velocity picks and prepare
% it for Tomography
clear; close all;

%fPathStore = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\FundVelFeb\AllFreqPicksPerc30\';
fPathStore = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\LowVelAug\SNR7\';

%load('grpTimesAllNew.mat');
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\FundVelFeb\GrpTimePickAllPerc30.mat')
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\LowVelAug\GrpTimePick.mat')
load('vLim.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% get the not-NaN indices from the travelTime matrix
fAll = 1.0:0.05:5;
snrCut = 7;
distCut = 600; % units in meters
lambdaCut = 1;
velMin = 0; velMax = 5000;
velEdges = velMin:20:velMax;
% update everything based on the stnDist
% for i = 1:1:length(rayAttribute(:,1))
%     if(rayAttribute(i,5)<=distCut)
%         grpTimeStore(:,i) = NaN;
%     end
% end

for fNo = 1:1:length(fAll)
    fSee = fAll(fNo);
    fStoreName = ['grpTimeNew',num2str(fSee),'.mat'];
    fSeeInd = find(fExt>=fSee,1,'first');
    
    notNaN = (~isnan(grpTimeStore(fNo,:)));
    %stnDistGood = (rayAttribute(:,5)>distCut);
    minDist = lambdaCut*vLim(fNo,2)/fSee;
    
    %%
    rayAttTemp = rayAttribute(notNaN,:);
    grpTimeStoreTemp = grpTimeStore(fNo,notNaN);
    grpTimeErrStoreTemp = grpTimeErrStore(fNo,notNaN);
    grpSNRStoreTemp = grpSNRStore(fSeeInd,notNaN);
    
    velDist = rayAttTemp(:,5)./grpTimeStoreTemp';
    
    % make the plots and get the stats
    %velMin = min(velDist)-100; velMax = max(velDist)+100;
    
    [velHist(:,fNo),edges] = histcounts(velDist,velEdges);
    velHist(velHist(:,fNo)==0,fNo) = NaN;
    velPrct50(fNo,1) = median(velDist);
    %%
    % now do the analysis after SNR cut-off
    goodSNRInd = (grpSNRStoreTemp>=snrCut);
    goodDistInd = (rayAttTemp(:,5)>=minDist);
    goodInd = goodSNRInd & goodDistInd';
    rayAttTemp = rayAttTemp(goodInd,:);
    grpTimeStoreTemp = grpTimeStoreTemp(goodInd);
    grpTimeErrStoreTemp = grpTimeErrStoreTemp(goodInd);
    
    velDistNew = rayAttTemp(:,5)./grpTimeStoreTemp';
    
    % now do the plots to check how everything works
    
    [velHistNew(:,fNo),~] = histcounts(velDistNew,velEdges);
    velHistNew(velHistNew(:,fNo)==0,fNo) = NaN;
    velPrct50New(fNo,1) = median(velDistNew);
    
    rayAttributeStoreFull = [];
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
end

edgesOut = (edges(1:end-1) + edges(2:end))/2;
figure(1);
subplot(1,2,1)
hold on;
surf(fAll,edgesOut,velHist);
shading interp; view(2);
plot3(fAll,velPrct50,10000*ones(length(fAll),1),'k','LineWidth',2);
hold off;

subplot(1,2,2)
hold on;
surf(fAll,edgesOut,velHistNew);
shading interp;view(2);
plot3(fAll,velPrct50New,10000*ones(length(fAll),1),'k','LineWidth',2);

hold off;