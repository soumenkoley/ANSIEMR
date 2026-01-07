% this script was written to pick the dispersion curve (phase velocity)
% March, 2022, S. Koley

clear; close all;
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\';
figHand = openfig([fPath,'subArray03\subArray03.fig']);
xlim([1,5]);
% now call getpts
[xi,yi] = getpts(figHand);