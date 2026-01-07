% this script was written to convert group to phse velocity
clear; %close all

fSt = 1.0;

% load the picked group velocity curve
%grpPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\GrpTimePickAllPerc30Median.mat';
%grpPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\StarsJan\OnlyFund\Star10\FundUpd.txt';
grpPath = 'A:\TestInver\StarsApril\Star50\FundUpdApril.txt';

%phPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\StarsJan\OnlyFund\Star1\PhFundAvg.txt';
phPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray00vFK.mat';

%loaded as vAll in workspace
load(phPath);
fStInd = find(vAll(:,1)>=fSt,1,'first');
% c0 = vAll(fStInd,2); % units in m/s
% k0 = 2*pi*fSt/c0;
c0 = 4000;
k0 = 2*pi*fSt/c0;

A = load(grpPath);
fStInd = find(A(:,1)==fSt,1,'first');
fAll = A(fStInd:end,1);
% fAll = fAll';
% integrate the group velocity with respect to frequency
k = 2*pi*cumtrapz(fAll(1:end),1./(1*A(fStInd:end,2)));
kC = k0-k(1,1);
k = kC+k;
phVelOut = 2*pi*fAll(1:end)./k;

% load the Dispersion curves from Terziet
% terzPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey2\assembleCrossCorr\RingSurveyDC\';
% obsDCBig = load([terzPath,'DCBigRing.txt']);
% obsDCSmall = load([terzPath,'DCSmallRing.txt']);
% 
% figure(1);
% hold on;
% plot3(fAll(:,1),phVelOut(:,1),1000*ones(length(fAll(:,1)),1),'ro','LineWidth',2);
% plot(obsDCBig(:,1),obsDCBig(:,2),'k','LineWidth',2);
% plot(obsDCSmall(:,1),obsDCSmall(:,2),'k','LineWidth',2);
% plot(vAll(:,1),vAll(:,2),'k','LineWidth',2);
% plot(vAll(:,1),vAll(:,2)+vAll(:,3),'k--','LineWidth',2);
% plot(vAll(:,1),vAll(:,2)-vAll(:,3),'k--','LineWidth',2);
% hold off;

% plot the results
figure(1011);
plot(vAll(:,1),vAll(:,2),'b','LineWidth',2);
hold on;
plot(fAll,phVelOut,'-o','color','r','LineWidth',2);
hold off;

figure(1012);
hold on;
plot3(fAll(:,1),phVelOut(:,1),1000*ones(length(fAll),1),'m','LineWidth',2);
hold off;