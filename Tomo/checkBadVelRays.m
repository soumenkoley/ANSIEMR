% this is just a test script
% I wanted to see which station pairs were giving apparent velocities
% that were much higher than the limits I had given

clear; close all;

% load teh final travel time picks
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\FundVelApril\GrpTimePickAllPerc20.mat');

% load the theoretical grpVals
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTimeFundApril.mat');
fShow = 1:0.05:5.0;

fStartInd = find(fExt>=fShow(1,1),1,'first');
fEndInd = find(fExt>=fShow(end),1,'first');

grpTimeStore = grpTimeStore(fStartInd:fEndInd,:);
for i = 1:1:length(rayAttribute(:,1))
    velGot(:,i) = rayAttribute(i,5)./grpTimeStore(:,i);
end


for i = 1:1:length(rayAttribute(:,1))
    badVelInd = [];
    stnDist = rayAttribute(i,5);
    velHigh = 1.2*stnDist./thTime(:,i);
    velLow = 0.8*stnDist./thTime(:,i);
    for j = 1:1:length(fShow)
        if((velGot(j,i)>velHigh(j,1)) || (velGot(j,i)<velLow(j,1)))
            badVelInd(j,1) = 1;
        end
        
    end
    if(~isempty(badVelInd))
        disp('Stop');
    end
end