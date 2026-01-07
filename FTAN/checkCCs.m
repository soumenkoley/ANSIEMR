% ok so I am writing this code to see all the cross correlation pairs in a
% subarray
% the idea is to select the ones with symmetric nature. Dont use others you
% will estimate the group velocity wong. Stick to the principles, dont be
% greedy :P

clear; close all;
dayVec1 = 10:1:30; dayVec2 = [];
dayVec = [dayVec1,dayVec2];
monthVec = [11*ones(1,length(dayVec1)),12*ones(1,length(dayVec2))];
yearVec = [20*ones(1,length(dayVec1)),20*ones(1,length(dayVec2))];

hrVec = [0:23];

dataPath = 'B:\LimburgBigSurvey1CC-Pair\';

%% SubArray S0
stnList = ["WUJ2A";"W3JWA";"W1KLA";"XBKUA";"XUKKA";"XDLDA";...
    "W1LKA";"X0K9A";"XHKAA";"W8KFA";"XUKWA";"XYKJA";"X2KQA";"YDK3A"]; % top left corner

%%
%stnList = ["0NQPA";"0ZQEA";"08QLA";"08P9A";"1MP1A";"1OP6A";"0KP6A";...
%    "0RPXA";"09PQA";"1PO5A";"0CPUA";"0PPFA";"0XPIA";"03PAA";"0UO6A";...
%    "0QO3A";"0GO1A";"0CO6A";"0POKA";"03OMA";"1EOGA"];

%stnList = ["X6NQA";"X7ORA";"Y0NZA";"Y3OLA";"YAOTA";"YEOMA";...
%    "YGNCA";"YHNNA";"YKNVA";"YLNBA";"YLOWA";"YONYA";"YRNEA";...
%    "YRNVA";"YSOJA";"YUNJA";"YUNZA";"YUOZA";"YXN9A";"YXOQA"];

% stnList = ["0CO6A";"0CPUA";"0GO1A";"Y3O5A";"Y7PAA";"Y8PQA";...
%     "Z2P2A";"Z3O2A";"Z5PHA";"Z6PPA";"ZCOUA";"ZGPDA";"ZGPVA";"ZIORA";...
%     "ZPPZA";"ZQO1A";"ZROJA";"ZRPCA";"ZXOSA";"ZYOMA";"ZZPAA"];

%stnList = ["ZROJA";"ZYOMA";"ZXOSA";"ZQO1A";"Z3O2A";"ZZPAA";"ZRPCA";...
%    "ZCOUA";"ZIORA";"ZGPDA"];

%stnList = ["0INXA";"0PNVA";"0XNZA";"0GNUA";"0GNDA";"0QNIA";"0ZNEA";...
%    "0PM7A";"0IM0A";"0POKA";"03OMA";"1EOGA";"1ANPA";"1BM8A";"1BMMA"];

%stnList = ["Z2MDA";"Z2NAA";"0PNVA";"0XNZA";"1ANPA";"1BM8A";"0ULZA";...
%    "0FL9A";"ZSMIA"];
%% S05 subarray
%stnList = ["Z2MDA";"0FL9A";"0ULZA";"0GK6A";"Z8LHA";"ZRLWA";"ZILIA";...
%    "Z2KWA";"Z6KQA";"ZUJ5A";"ZOKUA";"ZDK2A";"Y8KWA";"Y0KWA"];%S05

% SO5Full (+closely spaced sensors at Meleschet)
%stnList = ["Z2MDA";"0FL9A";"0ULZA";"0GK6A";"Z8LHA";"ZRLWA";"ZILIA";...
%    "Z2KWA";"Z6KQA";"ZUJ5A";"ZOKUA";"ZDK2A";"Y8KWA";"Y0KWA";"Y6MAA";...
%    "Y8MAA";"Y7L9B";"Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"]; 

%stnList = ["Z2KWA";"0GK6A";"Z6KQA";"ZUJ5A"];%S005
%stnList = ["Z2KWA";"ZOKUA";"ZDK2A";"ZILIA";"ZRLWA";"Z8LHA"]; %S0005
%% subArray S06
%stnList = ["YGNCA";"YRNEA";"YLNBA";"YUNJA";"YHNNA";"YKNVA";"YONYA";...
%    "YRNVA";"YUNZA";"Y0NZA";"X6NQA";"Y4NCA";"YGMXA";"X6M2A";...
%    "YLMPA";"YXN9A";"YSOJA";"YEOMA";"YLOWA";"YUOZA";"YXOQA";"Y3OLA";...
%    "ZLNPA";"ZCM4A";"YYMLA"];%S06

%% analyzing SO3 subarray better
%stnList = ["ZYOMA";"ZXOSA";"ZQO1A";"ZRPCA";"ZZPAA";"Z3O2A";"0CO6A";...
%    "Z5PHA";"0GO1A"];

%% subArray SO7, near left of Broek
%stnList = ["W5NUA";"W7N9A";"XDOGA";"XOOIA";"XTOIA";"XON9A";"XLN1A";...
%    "W5NUA";"XHNLA";"XTNZA";"XVNMA";"X6NQA";"XZNCA";"XHNAA";"XPM9A";...
%    "XVM7A";"XAOYA";"XNOYA";"XQO5A";"X3O0A";"YAOTA";"YEOMA"];
%% subArray S08
%stnList = ["YGKGA";"YHKVA";"YDK3A";"X2KQA";"X0K9A";"XXLUA";"X1LUA";...
%    "YHLWA";"X4MBA";"YEL9A";"XVMMA";"YGMXA";"YLMPA";"YYMLA";"ZBMKA";...
%    "Y0MBA";"Y8MAA";"ZILIA";"ZDK2A";"Y8KWA";"Y0KWA"];
%% subArray SO8A
%stnList = ["XUKKA";"XYKJA";"X2KQA";"XUKWA";"YGKGA";"YHKVA";"YDK3A";...
%    "X0K9A";"XXLUA";"X1LUA";"YHLWA";"YEL9A";"X4MBA";"YLMPA";"YGMXA";...
%    "X6M2A";"XMMBA";"XVMMA"];
%% subArray S08B
%stnList = ["Y0KWA";"Y8KWA";"ZDK2A";"ZILIA";"ZRLWA";"Z2MDA";"ZSMIA";...
%    "ZBMKA";"YYMLA";"Y0MBA";"Y6MAA";"Y7L9B";"Y8MAA";"Y7L9A";"Y8L8A";...
%    "Y9L7A";"Y8L6A";"YHLWA"];
%% S08-FULL
%stnList = ["X0K9A";"X1LUA";"X2KQA";"X4MBA";"X6M2A";"XMMBA";"XUKKA";"XUKWA";...
%    "XVMMA";"XXLUA";"XYKJA";"Y0KWA";"Y0MBA";"Y6MAA";"Y7L9A";"Y7L9B";...
%    "Y8KWA";"Y8L6A";"Y8L8A";"Y8MAA";"Y9L7A";"YDK3A";"YEL9A";"YGKGA";...
%    "YGMXA";"YHKVA";"YHLWA";"YLMPA";"YYMLA";"Z2MDA";"ZBMKA";"ZDK2A";...
%    "ZILIA";"ZRLWA";"ZSMIA"];
%% S09 is a new array, and includes the Meleschet closely spaced ones
 %stnList = ["Y0MBA";"YYMLA";"YLMPA";"YGMXA";"ZBMKA";"ZCM4A";...
 %    "Y4NCA";"Y0NZA";"YUNZA";"YONYA";"YKNVA";"YRNVA";"YHNNA";"YGNCA";...
 %    "YLNBA";"YRNEA";"YUNJA";"ZLNPA";"ZQNFA";"Y6MAA";"Y8MAA";"Y7L9B";...
 %    "Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"];
% this is S09 without Melleschet
%stnList = ["Y0MBA";"YYMLA";"YLMPA";"YGMXA";"ZBMKA";"ZCM4A";...
%    "Y4NCA";"Y0NZA";"YUNZA";"YONYA";"YKNVA";"YRNVA";"YHNNA";"YGNCA";...
%    "YLNBA";"YRNEA";"YUNJA";"ZLNPA";"ZQNFA";"Y8MAA"];
%% subarray 10, this is right of SO9, and includes the Melleschet
%stnList = ["Y0MBA";"YYMLA";"ZBMKA";"ZCM4A";"Y4NCA";"ZLNPA";"ZQNFA";...
%    "ZXNKA";"Z2NAA";"ZSMIA";"Z2MDA";"0FL9A";"0IM0A";"0GNDA";"ZRLWA";...
%    "Y6MAA";"Y8MAA";"Y7L9B";"Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A"];

% this is S10 but without Melleschet
%stnList = ["Y0MBA";"YYMLA";"ZBMKA";"ZCM4A";"Y4NCA";"ZLNPA";"ZQNFA";...
%    "ZXNKA";"Z2NAA";"ZSMIA";"Z2MDA";"0FL9A";"0IM0A";"0GNDA";"ZRLWA";...
%    "Y6MAA"];
% S10case1
%stnList = ["Y0MBA";"YYMLA";"ZBMKA";"ZCM4A";"Y4NCA";"ZLNPA";"ZQNFA";...
%    "ZXNKA";"Z2NAA";"ZSMIA";"Z2MDA";"0FL9A";"0IM0A";"0GNDA";"ZRLWA";...
%    "Y6MAA";"Y8MAA";"Y8L8A"];

%% subarray 11, ofcourse without Melleschet :P
%stnList = ["0ULZA";"0FL9A";"Z2MDA";"Z2NAA";"0GNDA";"0IM0A";"0PM7A";...
%    "0QNIA";"0ZNEA";"1BM8A";"1BMMA";"0GNUA";"0PNVA";"0XNZA";"0INXA"];

%% SubArray 12, overlaps with S11
%stnList = ["0IM0A";"Z2NAA";"ZXNKA";"0GNDA";"ZVN3A";"ZZN8A";"0POKA";...
%    "03OMA";"1EOGA";"1ANPA";"1BM8A";"0ZNEA";"0PM7A";"0QNIA";"0PNVA";...
%    "0XNZA";"0INXA";"0GNUA"];
%% subArray 13, overlaps with S0 and S09
%stnList = ["W4L7A";"XAL4A";"XMMBA";"XFMHA";"XELYA";"W9MRA";"XFMWA";...
%    "XVMMA";"XXLUA";"X1LUA";"X4MBA";"YEL9A";"YHLWA";"YLMPA";"YGMXA";...
%    "X6M2A";"XVM7A";"XPM9A";"XHNAA"]; % updated March 22, 2022

% this we call S13A
%stnList = ["W4L7A";"W9MRA";"XFMWA";"XFMHA";"XAL4A";...
%    "XELYA";"XMMBA";"XVMMA";"X4MBA";"XXLUA";"YHLWA";"YEL9A";...
%    "X1LUA";"XHNAA";"XPM9A";"XVM7A";"XZNCA";"X6M2A";"YGMXA";"YLMPA"];
%% SubArray 14, near Cottessen
% stnList = ["YCPBA";"X7PFA";"XYPWA";"X3QDA";"YIPMA";"YFPTA";"YAP0A";...
%     "YUP2A";"YOPFA";"X5QHA";"YAQGA";"YAQEA";"YCQGA";"YKQHA";"YIQEA";...
%     "YJQBA";"YHP8A";"YNP9A";"XPPNA";"XQO5A";"X3O0A";"YLOWA";"YUOZA";...
%     "Y8PQA";"Y7PAA"];
% subarray S14A
% stnList = ["YCPBA";"X7PFA";"XYPWA";"X3QDA";"YIPMA";"YFPTA";"YAP0A";...
%     "YUP2A";"YOPFA";"X5QHA";"YAQGA";"YAQEA";"YCQGA";"YKQHA";"YIQEA";...
%     "YJQBA";"YHP8A";"YNP9A"];
%% subArray 15, slightly northwest of S14, to see transition
%stnList = ["X3O0A";"X7PFA";"YCPBA";"YIPMA";"YOPFA";"YFPTA";"YUP2A";...
%    "Y8PQA";"Y7PAA";"Y3O5A";"Y3OLA";"YSOJA";"YXOQA";"YUOZA";"YLOWA";...
%    "YEOMA";"YAOTA";"X7ORA"];
%% subArray 16, left of S15
%stnList = ["XDOGA";"XAOYA";"XOOIA";"XTOIA";"XNOYA";"XQO5A";"X3O0A";...
%    "YCPBA";"X7PFA";"YIPMA";"YFPTA";"YAP0A";"XYPWA";"XPPNA";"XKPSA";...
%    "XIPOA";"XGPDA";"W7PSA"];

%% subArray 17, above S16 subArray
%stnList = ["W5NUA";"W7N9A";"XAOYA";"XDOGA";"XNOYA";"XQO5A";"X3O0A";...
%    "YAOTA";"X7ORA";"YEOMA";"X6NQA";"XVNMA";"XHNLA";"XLN1A";"XON9A";...
%    "XTOIA";"XOOIA";"XTNZA"];

% subArray17A, appepnding more stations to S17
%stnList = ["W5NUA";"W7N9A";"XAOYA";"XDOGA";"XNOYA";"XQO5A";"X3O0A";...
%    "YAOTA";"X7ORA";"YEOMA";"X6NQA";"XVNMA";"XHNLA";"XLN1A";"XON9A";...
%    "XTOIA";"XOOIA";"XTNZA";"YGNCA";"YHNNA";"YKNVA";"YONYA";"YUNZA";...
%    "Y0NZA";"YUNJA";"YRNEA";"YLNBA"];
%%
% trying small arrays now, to avoid heterogenity small1
%stnList = ["YGNCA";"YHNNA";"YKNVA";"YONYA";"YUNZA";"Y0NZA";"YRNVA";...
%    "YLNBA";"YRNEA";"YUNJA"];

% trying the array near Restaurand Buitenlust, small2
%stnList = ["X3O0A";"YCPBA";"YLOWA";"YEOMA";"YAOTA";"X7ORA";"YOPFA";"YUOZA";...
%    "YXOQA";"YSOJA"];

% a small array near Vakantie Vocken, that was super high velocity
%stnList = ["XON9A";"XOOIA";"XTOIA";"XTNZA";"XLN1A"];

% small array near Oud Bommerich
%stnList = ["XHNAA";"XPM9A";"XVM7A";"XZNCA";"XHNLA";"XVNMA";"XTNZA";...
%    "XLN1A"];

% another small array near HGN station
%stnList = ["W7PSA";"XIPOA";"XKPSA";"XGPDA";"XPPNA"];

% trying a small array within S13a, that was not very clear%
%stnList = ["X6M2A";"YGMXA";"YLMPA";"YEL9A";"X4MBA";"XVMMA";"XXLUA";...
%    "X1LUA";"YHLWA"];

% so now i try one with Melleschet, almost half of S13a
%stnList = ["Y6MAA";"Y8MAA";"Y7L9B";"Y7L9A";"Y8L8A";"Y9L7A";"Y8L6A";...
%    "Y0MBA";"YYMLA";"YLMPA";"YHLWA";"YEL9A";"YGMXA";"X1LUA";"X4MBA";...
%    "XXLUA";"XMMBA";"XVMMA";"X6M2A";"XUKKA";"XUKWA";"X0K9A";"XYKJA";...
%    "X2KQA";"YDK3A";"YHKVA";"YGKGA"];

% another small array within S0
%stnList = ["XUKKA";"XUKWA";"X0K9A";"XYKJA";"X2KQA";"YDK3A";...
%    "YHKVA";"YGKGA"]; % this is alright

% another small array left small6 array and within S13
%stnList= ["W4L7A";"XAL4A";"XMMBA";"XFMHA";"XELYA";"W9MRA";"XFMWA";...
%    "XVMMA";"XXLUA";"X1LUA";"X4MBA";"YEL9A";"YHLWA";"YLMPA";"YGMXA";...
%    "X6M2A";"XVM7A";"XPM9A";"XHNAA"];% small11 name

% another small array right of Epen
%stnList = ["W5NUA";"XHNLA";"XLN1A";"W7N9A";"XDOGA";"XON9A";"XTOIA";...
%    "XOOIA";"XNOYA";"XAOYA";"X6NQA";"X7ORA";"YEOMA";"X3O0A";"XGPDA";...
%    "XQO5A";"YAOTA"];

%stnList = ["ZOKUA";"ZRLWA"];
%stnList = ["XDLDA";"XBKUA"];
%stnList = ["XUKWA";"X0K9A";"X2KQA";"YDK3A"]; % S0, subArrayTest
%stnList = ["YGNCA";"YHNNA";"YKNVA";"YONYA";"YRNVA";"YUNJA";"YRNEA";"YLNBA"];% S1
%stnList = ["XMMBA";"XFMHA";"XVMMA";"W9MRA";"XFMWA";"XVM7A";...
%    "XPM9A";"XHNAA"]; % S2
%stnList = ["0PNVA";"0QNIA";"0XNZA";"1ANPA";"0ZNEA"];% S3
%stnList =["0GNUA";"0INXA";"0POKA";"03OMA";"0XNZA";"0PNVA"]; % S4
%stnList =["ZQNFA";"ZLNPA";"ZVN3A";"ZZN8A";"ZXNKA";"Z2NAA"];% S5
%stnList = ["Z2NAA";"ZXNKA";"ZVN3A";"ZZN8A";"0GNUA";"0INXA"];%S6
%stnList = ["ZXOSA";"ZQO1A";"ZRPCA";"ZZPAA";"Z3O2A";"ZYOMA";"ZROJA"];%S7
%stnList = ["ZZPAA";"Z5PHA";"Z6PPA";"0CPUA";"0GO1A";"0CO6A";"Z3O2A"];%S8
%stnList = ["YXOQA";"Y3O5A";"Y7PAA";"ZGPDA";"ZIORA";"ZCOUA";"Y3OLA";"YUOZA"];% S9
%stnList = ["YLOWA";"YUOZA";"YXOQA";"Y3OLA";"YXN9A";"YSOJA"];%S10
%stnList = ["0QNIA";"0PNVA";"0XNZA";"0ZNEA";"1ANPA"];
%stnList =["0XNZA";"03OMA";"0POKA";"0PNVA"];
%stnList = ["ZZPAA";"ZRPCA";"Z5PHA";"Z3O2A";"ZROJA";"ZQO1A";"0CO6A";...
%    "0GO1A";"ZXOSA";"ZYOMA";"0CPUA"];
%stnList = ["XVNMA";"X6NQA"];
%stnList = ["Y0MBA";"YHLWA";"Y6MAA";"Y7L9B";"Y7L9A";"Y8MAA";"Y8L8A";"Y9L7A";"Y8L6A"];
%stnList =  ["XYPWA";"YAP0A";"X3QDA";"X5QHA";"YAQGA";"YAQEA";"YCQGA";"YHP8A";"YJQBA";"YIQEA";"YKQHA"];
%stnList = ["YGNCA";"YLNBA";"YRNEA";"YUNJA";"YHNNA";"YKNVA";"YONYA";"YRNVA";"YUNZA";"Y0NZA"];
%stnList = ["XXLUA";"X1LUA";"YHLWA"];
%stnList = []

% Also read all the stn names
%stnList = unique(stnList);
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
stnLocs = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
refLong = 5.9060966; refLat = 50.7513927;

R = 6371000; % radius of earth in meters
for i = 1:1:length(stnList)
    stnInd = find(allStn==stnList(i,1));
    stnLocUse(i,1:2) = stnLocs(stnInd,1:2);
    stnLocUse(i,3) = stnInd;
    [xLong,yLat] = calculatedist(stnLocUse(i,1:2),[refLat,refLong],R);
    stnLocUse(i,4:5) = [xLong,yLat];
    %pairAz(i,3) = azimuth(pairNodeALatLong(1,1),pairNodeALatLong(1,2),pairNodeBLatLong(1,1),pairNodeBLatLong(1,2));
    
end

for i = 1:1:length(allStn)
    [xLong,yLat] = calculatedist(stnLocs(i,1:2),[refLat,refLong],R);
    allStnLoc(i,1:2) = [xLong,yLat]; 
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
for dayNo = 1:1:length(dayVec)
    if(dayVec(dayNo)<10)
        dayStr = ['0',num2str(dayVec(dayNo))];
    else
        dayStr = num2str(dayVec(dayNo));
    end
    
    foldPathA = [num2str(yearVec(dayNo)),num2str(monthVec(dayNo)),dayStr];
    
    for hrNo = 1:1:length(hrVec)
        ccStoreIndNow = 1;
        for stnNo = 1:1:length(stnList)
            stnA = stnList{stnNo,:};
            
            foldPathB = ['Hr',num2str(hrVec(hrNo))];
            totPath = [dataPath,'\',foldPathA,'\',foldPathB,'\'];
            fNameA = [stnA,'.mat'];
            
            fullPathA = [totPath,fNameA];
            if(exist(fullPathA))
                load(fullPathA);
                % convert character array to string array
                stnEnd = string(stnEnd);
                endStnList = setdiff(stnList,stnA,'stable');
                % loop for other stations
                for endStnNo = 1:1:length(endStnList)
                    stnInd = find(stnEnd==endStnList(endStnNo,1));
                    if(~isempty(stnInd))
                        ccFullStore(:,ccStoreIndNow) = ccFullStore(:,ccStoreIndNow) + ccStore(:,stnInd);
                        ccFullStoreNumber(ccStoreIndNow,1) = ccFullStoreNumber(ccStoreIndNow,1) + 1;
                    end
                %ccTemp = ccStore(:,stnInd);
                ccStoreIndNow = ccStoreIndNow+1;
                end
                
            else
                ccStoreIndNow = ccStoreIndNow+lList;
            end
            %disp('One stn done!')
        end
        %disp('One hour done!');
        fclose('all');
    end
    disp('One day done!');
end

% now time to sort the fullCCStore
ccCount = 1;
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
        ccStoreFinal(:,ccCount) = ccStoreFinal(:,ccCount)/ccStoreNumber(ccCount,1);
        if(isnan(ccStoreFinal(1,ccCount)))
            disp('STOP');
        end
        rayAttribute(ccCount,1:2) = stnLocUse(i,4:5);
        rayAttribute(ccCount,3:4) = stnLocUse(j,4:5);
        rayAttribute(ccCount,5) = sqrt((rayAttribute(ccCount,1)-rayAttribute(ccCount,3))^2+...
                       (rayAttribute(ccCount,2)-rayAttribute(ccCount,4))^2);
        rayAttribute(ccCount,6) = azimuth(stnLocUse(i,1),stnLocUse(i,2),stnLocUse(j,1),stnLocUse(j,2));
        rayAttribute(ccCount,7) = i; rayAttribute(ccCount,8) = j;
        ccCount = ccCount + 1;
    end
end

% remove NaN cross-corrs
ccStoreFinal = ccStoreFinal(:,(~isnan(ccStoreFinal(1,:))));
rayAttribute = rayAttribute((~isnan(ccStoreFinal(1,:))),:);
% now do one thing check the fft amplitudes and if they are all balanced
fftCCAll = fft(ccStoreFinal,1001,1);
fVec = linspace(0,1,1001)*fSamp;
figure(4)
plot(fVec,abs(fftCCAll));

figure(1)
%subplot(1,2,1)
plot(allStnLoc(:,1),allStnLoc(:,2),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;
plot(stnLocUse(:,4),stnLocUse(:,5),'ro','MarkerFaceColor','r','MarkerSize',6);

%stnSeeA = "W3JWA"; stnSeeB = "XHKAA";
stnSeeA = stnList(1,1); stnSeeB = stnList(2,1);

stnSeeAInd = find(stnList==stnSeeA);stnSeeBInd = find(stnList==stnSeeB);
plot([stnLocUse(stnSeeAInd,4);stnLocUse(stnSeeBInd,4)],[stnLocUse(stnSeeAInd,5);...
    stnLocUse(stnSeeBInd,5)],'k','LineWidth',2);
hold off;

stnASeeIndices = find(ccStoreNumber(:,2)==stnSeeAInd);
stnBSeeIndices = find(ccStoreNumber(:,3)==stnSeeBInd);
stnABIntersect = intersect(stnASeeIndices,stnBSeeIndices,'stable');

load('HdBp3_8Hz.mat');

% now see the cross-correlations one by one

goodCount = 1;
for i = 1:1:length(ccStoreFinal(1,:))
    ccNow = filtfilt(HdBp3_8Hz.Numerator,1,ccStoreFinal(:,i));
    plot(tArray,ccNow);
    prompt('Do you want to save this?')
    str = input(prompt,'s');
    
    if(str=='y')
        ccGood(:,goodCount) = ccStoreFinal(:,i);
        rayAttrGood(goodCount,:) = rayAttribute(i,:);
        goodCount = goodCount+1;
    end
end
