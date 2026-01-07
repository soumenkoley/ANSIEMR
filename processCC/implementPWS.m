% this script was written to implement PWS
% by S Koley, May 27, 2019

clear; close all;

fSamp = 25; % units in Hz
dt = 1/fSamp;

dayProc = 20:1:31; % month of Dec, 2017
monthNow = 'Dec';

% PWS parameters
K = 6;
nuStack = 2;

% CC Parameters
load('lagVal.mat');

% also load the rayAttribute file
load('E:\LimburgPassive2CCStore\RayAttribute.mat');

for ccNo = 83:1:83
    ccAllDays = [];
    tic;
    for dayNo = 1:1:length(dayProc)
        dayNow = dayProc(1,dayNo);
        if(dayNow<20)
            monthStr = 'Jan';
        else
            monthStr = 'Dec';
        end
        filePath = ['D:\LimburgPassiveHourlyCCs\',monthStr,num2str(dayNow),'\'];
        fileName = ['cc',num2str(ccNo),'Dec',num2str(dayNow),'.mat'];
        totFilePath = [filePath,fileName];
        load(totFilePath);
        %disp('file loaded');
        ccAllDays = [ccAllDays,ccPair];
    end
    % now get rid of the zeros
    [s1,s2] = size(ccAllDays);
    goodCCCount = 1;
    ccGoodHours = [];
    for allCCNo = 1:1:s2
        if(ccAllDays(1,allCCNo)~=0)
            ccGoodHours(:,goodCCCount) = ccAllDays(:,allCCNo);
            goodCCCount = goodCCCount +1;
        end
    end
    linearStack = mean(ccGoodHours,2);
    [PWSTrace,tVec] = PWS(ccGoodHours,linearStack,fSamp,K,nuStack,lagVal);
    
    figure(1)
    subplot(2,1,1)
    plot(tVec,linearStack);
    subplot(2,1,2)
    plot(tVec,PWSTrace);
    toc;
    disp('One CC done');
end