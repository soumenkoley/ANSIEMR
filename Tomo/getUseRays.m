% this code was written to get reliable raypaths
% we plot the ray paths with at least 2 Hz of group dispersion info

clear; close all;

load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\FundVelFeb\GrpTimePickAllPerc10.mat')
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');

% load the theoretical group times
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTimeFundFeb.mat');
nRays = length(rayAttribute(:,1));
fShow = 1.6:0.05:5;
fPlot = 1.6:0.05:5;
snrCut = 5;

fStartInd = find(fExt>=fShow(1),1,'first');
fEndInd = find(fExt>=fShow(end),1,'first');

grpTimeStore = grpTimeStore(fStartInd:fEndInd,:);
grpTimeErrStore = grpTimeErrStore(fStartInd:fEndInd,:);
grpSNRStore = grpSNRStore(fStartInd:fEndInd,:);


figure(1);
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ro','MarkerFaceColor','r',...
    'MarkerSize',6);
goodRayCount = 1;
for i = 128:1:nRays
    grpTimeStore(grpSNRStore(:,i)<=snrCut,i) = NaN;
    grpTimeErrStore(grpSNRStore(:,i)<=snrCut,i) = NaN;
    grpTimeStoreTemp = grpTimeStore(~isnan(grpTimeStore(:,i)),i);
    lGrp = length(grpTimeStoreTemp);
    if(lGrp>50)
%         figure(2);
%         plot(fShow,rayAttribute(i,5)./grpTimeStore(:,i));
%         hold on;
%         plot(fPlot,rayAttribute(i,5)./thTime(:,i),'r','LineWidth',2);
%         hold off;
        %disp('One check');
        figure(1);
        hold on;
        rayStartX = [rayAttribute(i,1);rayAttribute(i,3)];
        rayStartY = [rayAttribute(i,2);rayAttribute(i,4)];
        plot(rayStartX,rayStartY,'k','LineWidth',2);
        hold off;
        stnCheck(goodRayCount,1) = allStn(rayAttribute(i,7));
        stnCheck(goodRayCount,2) = allStn(rayAttribute(i,8));
        %disp([allStn(rayAttribute(i,7),1),' and ', allStn(rayAttribute(i,8),1)]);
        goodRayCount = goodRayCount+1;
    end
    
end
hold off;