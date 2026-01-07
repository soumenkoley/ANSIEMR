% this script was written to plot all the group velocity curves
% manually picked

clear; close all;

stnLocs = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
starPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\';

starNo = 1:1:81;

starExist = 1;
fExt = (1:0.05:5)';
meanVel = zeros(length(fExt),1);

figure(1);
hold on;
for i = 1:1:length(starNo)
    fPath = [starPath,'Star',num2str(i),'\FundUpdN.txt'];
    if(exist(fPath))
        load(fPath);
        
        plot(FundUpdN(:,1),FundUpdN(:,2),'k','LineWidth',2);
        fundInterp = interp1(FundUpdN(:,1),FundUpdN(:,2),fExt);
        meanVel = fundInterp+ meanVel;
        starLocPath = [starPath,'Star',num2str(i),'\starLocs.txt'];
        latLong = load(starLocPath);
        starLocs(starExist,1:2) = [latLong(1,1),latLong(2,1)];
        starExist = starExist+1;
    end
end
meanVel = meanVel/(starExist-1);
plot(FundUpdN(:,1),meanVel,'r','LineWidth',2);
hold off;

figure(2);
hold on;
plot(stnLocs(:,2),stnLocs(:,1),'bo','MarkerSize',6);
plot(starLocs(:,2),starLocs(:,1),'ko','MarkerSize',6,'MarkerFaceColor','k');
hold off;