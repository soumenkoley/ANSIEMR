% this script was written to implement PWS for one particular pair
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
load('RayAttribute.mat');

n1 = 156; n2 = 91;
n1Ind = find((rayAttribute(:,1)==n1) | (rayAttribute(:,1)==n2));
n2Ind = find((rayAttribute(:,2)==n1) | (rayAttribute(:,2)==n2));

pairInd = intersect(n1Ind,n2Ind);

if(~isempty(pairInd))
    nA = rayAttribute(pairInd,1);
    nB = rayAttribute(pairInd,2);
    
    if(nA == n1)
        % no need to flip
        flipInd = 0;
    else
        flipInd = 1;
    end
end

for ccNo = pairInd:1:pairInd
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
            
            if(flipInd)
                ccGoodHours(:,goodCount) = flipud(ccGoodHours(:,goodCount));
            end
            goodCCCount = goodCCCount +1;
        end
    end
    linearStack = mean(ccGoodHours,2);
    [PWSTrace,tVec] = PWS(ccGoodHours,linearStack,fSamp,K,nuStack,lagVal);
    
    figure(1)
    subplot(1,2,1)
    plot(tVec,linearStack);
    xlim([-10,10]);
    subplot(1,2,2)
    plot(tVec,PWSTrace);
    xlim([-10,10]);
    figure(2)
    plotseis(ccGoodHours,tVec,1:1:(goodCCCount-1));
    ylim([-6,6]);
    toc;
    disp('One CC done');
end