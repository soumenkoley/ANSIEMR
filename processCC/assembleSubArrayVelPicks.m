% this script was written to asemble all the phase velocity pics from all
% subArrays
clear; close all;
fPath = 'B:\LimburgBigSurvey1CC-Pair\SubArrayVelPicks\';

arrayNum = 0:1:20;

rayAttributeStoreFull = [];
phaseVelStoreFull = [];
timePickStoreFull = [];
errorTStoreFull = [];
for i = 1:1:length(arrayNum)
    if(arrayNum(i)<10)
        arrayStr = ['0',num2str(arrayNum(i))];
    else
        arrayStr = num2str(arrayNum(i));
    end
    fPathTot = [fPath,'subArray',arrayStr,'Grp.mat'];
    if(exist(fPathTot))
        load(fPathTot);
        
        indGood = find(rayAttributeStore(3,:)>300);
        rayAttributeStore = rayAttributeStore(:,indGood);
        phaseVelStore = phaseVelStore(:,indGood);
        timePickStore = timePickStore(:,indGood);
        errorTStore = errorTStore(:,indGood);
        
        rayAttributeStoreFull = [rayAttributeStoreFull,rayAttributeStore];
        timePickStoreFull = [timePickStoreFull,timePickStore];
        phaseVelStoreFull = [phaseVelStoreFull,phaseVelStore];
        errorTStoreFull = [errorTStoreFull,errorTStore];
        disp('file loaded');
    end
end

% eliminate redundant ray-paths
for i = 1:1:length(rayAttributeStoreFull(1,:))
    rayId(i) = min(rayAttributeStoreFull(1,i),rayAttributeStoreFull(2,i))*1000+...
        max(rayAttributeStoreFull(1,i),rayAttributeStoreFull(2,i));
end
[a,uniqueInd] = unique(rayId,'stable');
rayAttributeStoreFull = rayAttributeStoreFull(:,uniqueInd);
phaseVelStoreFull = phaseVelStoreFull(:,uniqueInd);
timePickStoreFull = timePickStoreFull(:,uniqueInd);
errorTStoreFull = errorTStoreFull(:,uniqueInd);