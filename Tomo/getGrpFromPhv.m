% this script was written to get group velocity from the beamforming
% phase velocity picks
clear; %close all
% load the fundamental and the overtone
%fund = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\FundNew.txt');

%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray09vFKNew.mat');
vCCBF = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\modeBFK.txt');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\vFKHighPickFullArray.mat');

figure(184);
plot(vCCBF(:,1),vCCBF(:,2),'b','LineWidth',2);
hold on;
plot(vAll(:,1),vAll(:,2),'r','LineWidth',2);
hold off;
legend({'CCBF','FK'});
%vAll = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\StarsJan\OnlyFund\Star1\PhFundAvg - Copy.txt');

fundCCBF = vCCBF;
fundCCBF(:,2) = smooth(fundCCBF(:,2));

fundFK = vAll;
fundFK(:,2) = smooth(fundFK(:,2));

[fundGrpCCBF] = convPh2Grp(fundCCBF(:,1),fundCCBF(:,2));
[fundGrpFK] = convPh2Grp(fundFK(:,1),fundFK(:,2));

%overt = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\OvertNew.txt');

% [overtGrp] = convPh2Grp(overt(:,1),overt(:,2));

% figure(182);
% hold on;
% plot(fundCCBF(:,1),fundCCBF(:,2),'b','LineWidth',2);
% plot(fundFK(:,1),fundFK(:,2),'r','LineWidth',2);
% hold off;

figure(183);
hold on;
plot(fundGrpCCBF(:,1),fundGrpCCBF(:,2),'b','LineWidth',2);
plot(fundGrpFK(:,1),fundGrpFK(:,2),'r','LineWidth',2);
%plot(overtGrp(:,1),overtGrp(:,2),'r','LineWidth',2);
hold off;
legend({'CCBF','FK'});