% this script was writte to assemble the cross-correlations generated from
% the python code getCrossCorr.ipynb
% we generate only cross-corr pairs

clear; %close all;
subPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\';
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');

%dataPath = 'B:\LimburgBigSurvey1CC-Pair\';
dataPath = 'D:\LimburgSurvey1CC\stack\';
%dataPath = 'D:\ETSurvey1a\';

%%
% load the overtone or the fundamental derived in that subArray
velUse = load([subPath,'subArray20\OvertGrp.txt']);
storeStr = 'subArray20Grp';
fExt = 1.2:0.1:2.5;
fExtIndA = find(velUse(:,1)>=fExt(1,1),1,'first');
fExtIndB = find(velUse(:,1)>=fExt(1,end),1,'first');
velUse = velUse(fExtIndA:fExtIndB,:);

% error in observations
minErr = -1;
maxErr = 1; % in percentage

%% SubArray S0
% stnList = ["WUJ2A";"W3JWA";"W1KLA";"XBKUA";"XUKKA";"XDLDA";...
%    "W1LKA";"X0K9A";"XHKAA";"W8KFA";"XUKWA";"XYKJA";"X2KQA";...
%    "YDK3A";"YGKGA";"YHKVA";"W4L7A";"XAL4A";"XELYA";"XXLUA";"X1LUA"];
% improved low-freq resolution
%% subArray S01
% stnList = ["0NQPA";"0ZQEA";"08QLA";"08P9A";"1MP1A";"1OP6A";"0KP6A";...
%    "0RPXA";"09PQA";"1PO5A";"0CPUA";"0PPFA";"0XPIA";"03PAA";"0UO6A";...
%    "0QO3A";"0GO1A";"0CO6A";"0POKA";"03OMA";"1EOGA"];

%% subArray S02, on top of S03
% stnList = ["ZCM4A";"Y4NCA";"YUNJA";"YRNEA";"YRNVA";"YUNZA";"Y0NZA";...
%    "YXN9A";"YSOJA";"YXOQA";"Y3OLA";"ZCOUA";"ZIORA";"ZQO1A";"ZXOSA";...
%    "ZYOMA";"ZROJA";"ZZN8A";"ZVN3A";"ZXNKA";"ZLNPA";"ZQNFA"];

%% S05 subarray
% SO5Full (+closely spaced sensors at Meleschet)
% stnList = ["Z2MDA";"0FL9A";"0ULZA";"0GK6A";"Z8LHA";"ZRLWA";"ZILIA";...
%    "Z2KWA";"Z6KQA";"ZUJ5A";"ZOKUA";"ZDK2A";"Y8KWA";"Y0KWA";"Y6MAA";...
%    "Y8MAA";"Y7L9B";"Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"];

%% analyzing SO3 subarray better
% stnList = ["0CO6A";"0CPUA";"0GO1A";"Y3O5A";"Y7PAA";"Y8PQA";"Z2P2A";...
%    "Z3O2A";"Z5PHA";"Z6PPA";"ZCOUA";"ZGPDA";"ZGPVA";"ZIORA";"ZPPZA";...
%    "ZQO1A";"ZROJA";"ZRPCA";"ZXOSA";"ZYOMA";"ZZPAA"];

%% S08-FULL
% stnList = ["X0K9A";"X1LUA";"X2KQA";"X4MBA";"X6M2A";"XMMBA";"XUKKA";"XUKWA";...
%    "XVMMA";"XXLUA";"XYKJA";"Y0KWA";"Y0MBA";"Y6MAA";"Y7L9A";"Y7L9B";...
%    "Y8KWA";"Y8L6A";"Y8L8A";"Y8MAA";"Y9L7A";"YDK3A";"YEL9A";"YGKGA";...
%        "YGMXA";"YHKVA";"YHLWA";"YLMPA";"YYMLA";"Z2MDA";"ZBMKA";"ZDK2A";...
%    "ZILIA";"ZRLWA";"ZSMIA"];
%% S09 is a new array, and includes the Meleschet closely spaced ones
% stnList = ["Y0MBA";"YYMLA";"YLMPA";"YGMXA";"ZBMKA";"ZCM4A";...
%     "Y4NCA";"Y0NZA";"YUNZA";"YONYA";"YKNVA";"YRNVA";"YHNNA";"YGNCA";...
%     "YLNBA";"YRNEA";"YUNJA";"ZLNPA";"ZQNFA";"Y6MAA";"Y8MAA";"Y7L9B";...
%     "Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"];
%% subarray 10, this is right of SO9, and includes the Melleschet
% stnList = ["Y0MBA";"YYMLA";"ZBMKA";"ZCM4A";"Y4NCA";"ZLNPA";"ZQNFA";...
%    "ZXNKA";"Z2NAA";"ZSMIA";"Z2MDA";"0FL9A";"0IM0A";"0GNDA";"ZRLWA";...
%    "Y6MAA";"Y8MAA";"Y7L9B";"Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"];

%% subarray 11, ofcourse without Melleschet :P
% stnList = ["0ULZA";"0FL9A";"Z2MDA";"Z2NAA";"0GNDA";"0IM0A";"0PM7A";...
%    "0QNIA";"0ZNEA";"1BM8A";"1BMMA";"0GNUA";"0PNVA";"0XNZA";"0INXA"];

%% SubArray 12, overlaps with S11
% stnList = ["0IM0A";"Z2NAA";"ZXNKA";"0GNDA";"ZVN3A";"ZZN8A";"0POKA";...
%    "03OMA";"1EOGA";"1ANPA";"1BM8A";"0ZNEA";"0PM7A";"0QNIA";"0PNVA";...
%    "0XNZA";"0INXA";"0GNUA"];
%% subArray 13, overlaps with S0 and S09
% stnList = ["W4L7A";"XAL4A";"XELYA";"W9MRA";"XFMHA";"XMMBA";"XFMWA";...
%    "XVMMA";"XHNAA";"XPM9A";"XVM7A"];
%% SubArray 14, near Cottessen
% stnList = ["YCPBA";"X7PFA";"XYPWA";"X3QDA";"YIPMA";"YFPTA";"YAP0A";...
%     "YUP2A";"YOPFA";"X5QHA";"YAQGA";"YAQEA";"YCQGA";"YKQHA";"YIQEA";...
%     "YJQBA";"YHP8A";"YNP9A";"XPPNA";"X3O0A";"YAOTA";"X7ORA";"YEOMA";...
%     "YLOWA";"YUOZA"];

%% subArray 15, slightly northwest of S14, to see transition
% stnList = ["X3O0A";"X7PFA";"YCPBA";"YIPMA";"YOPFA";"YFPTA";"YUP2A";...
%    "Y8PQA";"Y7PAA";"Y3O5A";"Y3OLA";"YSOJA";"YXOQA";"YUOZA";"YLOWA";...
%    "YEOMA";"YAOTA";"X7ORA"];
%% subArray 17, above S16 subArray
% stnList = ["W5NUA";"W7N9A";"XAOYA";"XDOGA";"XNOYA";"XQO5A";"X3O0A";...
%    "YAOTA";"X7ORA";"YEOMA";"X6NQA";"XVNMA";"XHNLA";"XLN1A";"XON9A";...
%    "XTOIA";"XOOIA";"XTNZA"];
%% subArray S18, above S15, overlaps with S09 annd S15
% stnList = ["YUNJA";"YRNEA";"YLNBA";"YGNCA";"YHNNA";"YKNVA";"YONYA";...
%   "YRNVA";"YUNZA";"Y0NZA";"YXN9A";"Y3OLA";"YXOQA";"YSOJA";"YLOWA";...
%   "YEOMA";"YUOZA";"X6NQA";"YAOTA";"X7ORA";"X6M2A";"XZNCA";"XVNMA";...
%   "XTNZA";"XOOIA";"XTOIA";"XON9A"];
%% subArray 19, near Terziet
% stnList = ["WKQQA";"W3QZA";"W2RAA";"XFRFA";"W4QGA";"XPQJA";"W7PSA";...
%    "XKPSA";"XGPDA";"XIPOA";"XPPNA";"XYPWA";"X3QDA";"X5QHA";"YCPBA";...
%    "X7PFA";"XYPWA";"X3QDA";"YIPMA";"YFPTA";"YAP0A";"YUP2A";"YOPFA";...
%    "X5QHA";"YAQGA";"YAQEA";"YCQGA";"YKQHA";"YIQEA";"YJQBA";"YHP8A";...
%    "YNP9A"];
%% subArray 20, this bridges gap between S12, S02, S03
stnList = ["ZQNFA";"ZXNKA";"Z2NAA";"0GNDA";"0QNIA";"0PNVA";"0INXA";...
    "0GNUA";"ZZN8A";"ZVN3A";"ZROJA";"ZYOMA";"ZXOSA";"Z3O2A";"0CO6A";...
    "0GO1A";"0POKA";"0XNZA";"03OMA";"ZLNPA"];
%%
% Also read all the stn names
%stnList = unique(stnList);
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
%stnLocs = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
stnList = unique(stnList);
for i = 1:1:length(stnList)
    stnInd = find(allStn(:,1)==stnList(i,1));
    stnLocs(i,1) = stnInd;
    stnLocs(i,2:3) = nodeLocationsCartesian(stnInd,2:3);
end


% create an array of indices for storing the CCs
lList = length(stnList)-1;
for i = 1:1:(lList+1)
    ccStoreInd(((i-1)*lList+1):(i*lList),1) = i; 
    ccStoreInd(((i-1)*lList+1):(i*lList),2) = [(1:(i-1))';((i+1):(lList+1))'];
end
% Also read all the stn names

nCC = (lList+1)*lList/2;
ccFullStore = zeros(1001,((lList+1)*(lList)));
ccFullStoreNumber = zeros(((lList+1)*(lList)),1);

fSamp = 25;
tArray = [(-500:1:-1),(0:1:500)]/fSamp;
% now loop
ccStoreIndNow = 1;
for stnNo = 1:1:length(stnList)
    stnA = stnList{stnNo,:};
    
    fNameA = [stnA,'.mat'];
    
    fullPathA = [dataPath,fNameA];
    if(exist(fullPathA))
        load(fullPathA);
        % convert character array to string array
        stnEnd = string(stnStore);
        endStnList = setdiff(stnList,stnA,'stable');
        % loop for other stations
        for endStnNo = 1:1:length(endStnList)
            stnInd = find(stnEnd==endStnList(endStnNo,1));
            if(~isempty(stnInd))
                ccFullStore(:,ccStoreIndNow) = ccFullStore(:,ccStoreIndNow) + ccStoreFinal(:,stnInd);
            end
            %ccTemp = ccStore(:,stnInd);
            ccStoreIndNow = ccStoreIndNow+1;
        end
        
    else
        ccStoreIndNow = ccStoreIndNow+lList;
    end
    %disp('One stn done!')
end

% now time to sort the fullCCStore
ccCount = 1;
ccStoreFinal = [];
stnStore = strings(nCC,2);
for i = 1:1:lList
    for j = (i+1):1:(lList+1)
        % case 1
        
        indA = find(ccStoreInd(:,1)==i);
        indB = find(ccStoreInd(:,2)==j);
        indAB = intersect(indA,indB,'stable');
        ccStoreFinal(:,ccCount) = ccFullStore(:,indAB);
        ccStoreNumber(ccCount,1) = ccFullStoreNumber(indAB,1);
        ccStoreNumber(ccCount,2:3) = [i,j];
        indA = find(ccStoreInd(:,2)==i);
        indB = find(ccStoreInd(:,1)==j);
        indBA = intersect(indA,indB,'stable');
        ccStoreFinal(:,ccCount) = ccStoreFinal(:,ccCount) + flipud(ccFullStore(:,indBA));
        ccStoreNumber(ccCount,1) = ccStoreNumber(ccCount,1) + ccFullStoreNumber(indBA,1);
        
        if(isnan(ccStoreFinal(1,ccCount)))
            disp('STOP');
        end
        if(~isempty(indAB))
            rayAttributeStore(1,ccCount) = stnLocs(ccStoreInd(indAB,1),1);
            rayAttributeStore(2,ccCount) = stnLocs(ccStoreInd(indAB,2),1);
            stn1 = stnLocs(ccStoreInd(indAB,1),2:3);
            stn2 = stnLocs(ccStoreInd(indAB,2),2:3);
            rayAttributeStore(3,ccCount) = sqrt((stn1(1,1)-stn2(1,1))^2+...
                (stn1(1,2)-stn2(1,2))^2);
            rayAttributeStore(4:5,ccCount) = [stn1(1,1);stn1(1,2)];
            rayAttributeStore(6:7,ccCount) = [stn2(1,1);stn2(1,2)];
            rayAttributeStore(8,ccCount) = findSlopeNew(stn1,stn2);
        end
        
        if(~isempty(indBA))
            rayAttributeStore(1,ccCount) = stnLocs(ccStoreInd(indBA,1),1);
            rayAttributeStore(2,ccCount) = stnLocs(ccStoreInd(indBA,2),1);
            stn1 = stnLocs(ccStoreInd(indBA,1),2:3);
            stn2 = stnLocs(ccStoreInd(indBA,2),2:3);
            rayAttributeStore(3,ccCount) = sqrt((stn1(1,1)-stn2(1,1))^2+...
                (stn1(1,2)-stn2(1,2))^2);
            rayAttributeStore(4:5,ccCount) = [stn1(1,1);stn1(1,2)];
            rayAttributeStore(6:7,ccCount) = [stn2(1,1);stn2(1,2)];
            rayAttributeStore(8,ccCount) = findSlopeNew(stn1,stn2);
        end
        stnStore(ccCount,1) = allStn(stnLocs(i,1),1);
        stnStore(ccCount,2) = allStn(stnLocs(j,1),1);
        rng(ccCount,'twister');
        errVec = (maxErr-minErr).*rand(length(fExt),1) + minErr;
        phaseVelStore(:,ccCount) = velUse(:,2) + velUse(:,2).*(errVec/100);
        timePickStore(:,ccCount) = (rayAttributeStore(3,ccCount)./velUse(:,2));
        errorTStore(:,ccCount) = timePickStore(:,ccCount).*(errVec/100);
        ccCount = ccCount + 1;
    end
end
% remove NaN cross-corrs
goodInd = ~isnan(ccStoreFinal(1,:));
ccStoreFinal = ccStoreFinal(:,goodInd);
rayAttributeStore = rayAttributeStore(:,goodInd);
phaseVelStore = phaseVelStore(:,goodInd);
timePickStore = timePickStore(:,goodInd);
errorTStore = errorTStore(:,goodInd);
stnStore = stnStore(goodInd,:);

figure(1)
hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo',...
    'MarkerFaceColor','b','MarkerSize',8);
plot(stnLocs(:,2),stnLocs(:,3),'ro',...
    'MarkerFaceColor','k','MarkerSize',8);

storePath = ['B:\LimburgBigSurvey1CC-Pair\SubArrayVelPicks\',storeStr,'.mat'];

save(storePath,'ccStoreFinal','errorTStore','fExt','phaseVelStore',...
    'rayAttributeStore','stnList','stnLocs','tArray','timePickStore');

