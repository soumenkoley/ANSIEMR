% this scrpt uses the phase velocity curves derived from FK analysis
% to get the mean phase velocity and the upper and lower limits
% to be used for inversion

clear; close all;

vFKPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\';

load([vFKPath,'subArray00vFK.mat']);
vAssem(:,1) = vAll(:,2);
load([vFKPath,'subArray04vFK.mat']);
vAssem(:,2) = vAll(:,2);
load([vFKPath,'subArray05vFK.mat']);
vAssem(:,3) = vAll(:,2);
load([vFKPath,'subArray0511vFK.mat']);
vAssem(:,4) = vAll(:,2);
% load([vFKPath,'subArray08CvFK.mat']);
% vAssem(:,5) = vAll(:,2);
fVec = vAll(:,1);
vMean(:,1) = fVec;
% get the mean phase velocity
vMean(:,2) = mean(vAssem,2);


% get the upper and lower limits
for i = 1:1:length(fVec)
    vLimits(i,1) = min(vAssem(i,:));
    vLimits(i,2) = max(vAssem(i,:));
end

% get the max of difference between vMean and vLimits

for i = 1:1:length(fVec)
    v1 = abs(vMean(i,2) - vLimits(i,1));
    v2 = abs(vMean(i,2) - vLimits(i,2));
    vErr(i,1) = max(v1,v2);
    
end

vMean(:,3) = vErr;

figure(1); subplot(1,2,1); hold on;
plot(vAll(:,1),vAssem(:,1),'b','LineWidth',2);
plot(vAll(:,1),vAssem(:,2),'r','LineWidth',2);
plot(vAll(:,1),vAssem(:,3),'g','LineWidth',2);
plot(vAll(:,1),vAssem(:,4),'m','LineWidth',2);
%plot(vAll(:,1),vAssem(:,5),'c','LineWidth',2);
legend({'subArray 00','SubArray 04','SubArray 05','SubArray 11' });
hold off;

figure(1); subplot(1,2,2); hold on;
plot(vAll(:,1),vMean(:,2),'k','LineWidth',2);
plot(vAll(:,1),vLimits(:,1),'k--','LineWidth',2);
plot(vAll(:,1),vLimits(:,2),'k--','LineWidth',2);

% load the fundamental from FK analysis of all arrays 09, 06, 03, 02, 15,
% 08C
%load([vFKPath,'vFKFundPickAllArrays.mat']);
load([vFKPath,'subArray00040511vFK.mat'])
plot(vAll(:,1),vAll(:,2),'r','LineWidth',2);
hold off;