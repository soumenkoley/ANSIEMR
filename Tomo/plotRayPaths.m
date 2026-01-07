% this script was written to load a subArray phase velocity file and plot
% the ray paths
clear; %close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
fPath = 'B:\LimburgBigSurvey1CC-Pair\SubArrayVelPicks\';
subArrayName = 'subArray02.mat';
% load all station location file
load('nodeLocationsCartesian.mat');

load([fPath,subArrayName]);
figure(1);
hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo',...
    'MarkerFaceColor','b');
for i = 1:1:length(rayAttributeStore(1,:))
     plot([rayAttributeStore(4,i),rayAttributeStore(6,i)],[rayAttributeStore(5,i),rayAttributeStore(7,i)],'k-');
end
hold off;