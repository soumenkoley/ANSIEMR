% this script was written to see the CCs at a fixed separation, but
% changing azimuth
clear; close all;

fPath = 'B:\LimburgBigSurvey1CC-Pair\SubArrayCC\';

% load all the CCs,
load([fPath,'allCCs.mat']);

stnSep = 2500;
stnSepErr = 100;
stnDist1 = stnSep-stnSepErr;
stnDist2 = stnSep+stnSepErr;

% find the indexes of the stations from ray attribute
indStn1 = find(rayAttribute(:,5)>=stnDist1);
indStn2 = find(rayAttribute(:,5)<=stnDist2);
indStn = intersect(indStn1,indStn2);

ccStore = ccStoreFinal(:,indStn);
ccAttribute = rayAttribute(indStn,:);

load('HdBp1_4Hz.mat');
for i = 1:1:length(indStn)
    ccStore(:,i) = filtfilt(HdBp1_4Hz.Numerator,1,ccStore(:,i));
    %ccStore(:,i) = ccStore(:,i)/max(abs(ccStore(:,i)));
end

% add the CCs in 10^0 azimuths
dAz = 8;
azBins = 0:dAz:360;

for i = 1:1:(length(azBins)-1)
    azStart = azBins(i); azEnd = azBins(i+1);
    azInd1 = find(ccAttribute(:,6)>=azStart);
    azInd2 = find(ccAttribute(:,6)<=azEnd);
    azInd = intersect(azInd1,azInd2);
    ccStoreNew(:,i) = mean(ccStore(:,azInd),2);
    ccStoreNew(:,i) = ccStoreNew(:,i)/max(abs(ccStoreNew(:,i)));
end
figure(1)
imagesc(azBins(1:(end-1)),tArray,ccStoreNew);
colorbar;colormap('jet');
ylim([-6,6])

figure(2)
plotseis(ccStoreNew,tArray,azBins(1:(end-1)));
camroll(90);