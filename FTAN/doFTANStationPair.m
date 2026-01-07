% doFTAN
clear;
% this function performs FTAN on a desired station pair
% just wrote it for myself for analyzing how the data looks :)

% first load the rayAttribute file
addpath('C:\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\')
load('D:\GreensFunction\psdCorrAvgDay4-15.mat');
load('D:\GreensFunction\rayAttributeDay4-15.mat');
load('D:\GreensFunction\lagValDay4-15.mat');
n1 = 177; n2 = 94;
fSamp = 25; % remember you downsampled it from 250-25 Hz
fExt = (2:0.1:8)';
plotSet = 1;

n1Ind = find((rayAttribute(:,1)==n1) | (rayAttribute(:,1)==n2));
n2Ind = find((rayAttribute(:,2)==n1) | (rayAttribute(:,2)==n2));

pairInd = intersect(n1Ind,n2Ind);

if(~isempty(pairInd))
    nA = rayAttribute(pairInd,1);
    nB = rayAttribute(pairInd,2);
    
    if(nA == n1)
        % no need to flip
        ccPair = psdCorrAvg(:,pairInd);
    else
        ccPair = -flipud(psdCorrAvg(:,pairInd));
    end
end

[SAIFFT,fExt,distOut,timePickUpdated,velPickUpdated,timePickResolve] = ...
                                          FTAN(ccPair,lagVal,200,fSamp,fExt,n1,n2,plotSet);