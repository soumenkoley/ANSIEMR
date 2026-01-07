% this is an implementation of phase weighted stacking following Scheimmel
% et al 2011
% S Koley, March, 2022
clear; close all;
dayVec1 = [13:1:20]; dayVec2 = [1:2];
dayVec = [dayVec1,dayVec2];
monthVec = [11*ones(1,length(dayVec1)),12*ones(1,length(dayVec2))];
yearVec = [20*ones(1,length(dayVec1)),20*ones(1,length(dayVec2))];

hrVec = 0:1:23;

dataPath = 'B:\LimburgBigSurvey1CC-Pair\';

%stnA = 'ZCM4A'; stnB = 'ZQNFA';
%stnList = ["YRNVA";"YHNNA"];
stnList = ["YGNCA";"YHNNA";"YKNVA";"YONYA";"YUNZA";"Y0NZA";"YRNVA";...
    "YLNBA";"YRNEA";"YUNJA"];
%stnList = ["Y0MBA";"YYMLA";"ZBMKA";"ZCM4A";"Y4NCA";"ZLNPA";"ZQNFA";...
%    "ZXNKA";"Z2NAA";"ZSMIA";"Z2MDA";"0FL9A";"0IM0A";"0GNDA";"ZRLWA";...
%    "Y6MAA"];

% Also read all the stan names
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

ccAvg = zeros(1001,1);
ccFullStore = [];
fSamp = 25;
% try phase weighted stacking
K = 6; nuStack = 3;
lagVal = (-500:1:500);
tArray = lagVal/fSamp;
% now loop
ccCount = 1;

load('HdBp3_8Hz.mat');
ccStack = 1;
for stnANo = 1:1:(length(stnList)-1)
    for stnBNo = (stnANo+1):1:length(stnList)
        ccCount = 1;
        stnA = stnList(stnANo,1); stnB = stnList(stnBNo,1);
        % convert to char
        stnA = convertStringsToChars(stnA);
        stnB = convertStringsToChars(stnB);
        for dayNo = 1:1:length(dayVec)
            if(dayVec(dayNo)<10)
                dayStr = ['0',num2str(dayVec(dayNo))];
            else
                dayStr = num2str(dayVec(dayNo));
            end
            
            foldPathA = [num2str(yearVec(dayNo)),num2str(monthVec(dayNo)),dayStr];
            
            for hrNo = 1:1:length(hrVec)
                foldPathB = ['Hr',num2str(hrVec(hrNo))];
                totPath = [dataPath,'\',foldPathA,'\',foldPathB,'\'];
                fNameA = [stnA,'.mat'];
                
                fullPathA = [totPath,fNameA];
                if(exist(fullPathA))
                    load(fullPathA);
                    stnInd = getStnInd(stnEnd,stnB);
                    %ccTemp = ccStore(:,stnInd);
                    if(stnInd==0)
                        fNameB = [stnB,'.mat'];
                        fullPathB = [totPath,fNameB];
                        if(exist(fullPathB))
                            load(fullPathB);
                            stnInd = getStnInd(stnEnd,stnA);
                            if(stnInd==0)
                                ccTemp = zeros(1001,1);
                            else
                                ccTemp = flipud(ccStore(:,stnInd));
                            end
                        else
                            ccTemp = zeros(1001,1);
                        end
                    else
                        ccTemp = ccStore(:,stnInd);
                    end
                    
                    if(ccTemp(1,1)~=0)
                        ccFullStore(:,ccCount) = ccTemp;
                        ccCount = ccCount+1;
                    end
                else
                    disp('StnA not found');
                    break;
                end
            end
        end
        % now also store the rayAttribute
        rayAttribute(ccStack,1:2) = stnLocUse(stnANo,4:5);
        rayAttribute(ccStack,3:4) = stnLocUse(stnBNo,4:5);
        rayAttribute(ccStack,5) = sqrt((rayAttribute(ccStack,1)-rayAttribute(ccStack,3))^2+...
                       (rayAttribute(ccStack,2)-rayAttribute(ccStack,4))^2);
        rayAttribute(ccStack,6) = azimuth(stnLocUse(stnANo,1),stnLocUse(stnANo,2),...
            stnLocUse(stnBNo,1),stnLocUse(stnBNo,2));
        rayAttribute(ccStack,7) = stnANo; rayAttribute(ccStack,8) = stnBNo;
        disp('One pair read!');
        % Now do PWS
        for i = 1:1:length(ccFullStore(1,:))
            %ccFullStore(:,i) = filtfilt(HdBp3_8Hz.Numerator,1,ccFullStore(:,i));
        end
        linearStack(:,ccStack) = mean(ccFullStore,2);
        analyticSig = hilbert(ccFullStore);
        phaseNow = angle(analyticSig);
        sqrtUnity = sqrt(-1);
        phaseStack = mean(exp(sqrtUnity*phaseNow),2);
        phaseStack = (abs(phaseStack)).^2;
        phaseMat = repmat(phaseStack,1,length(ccFullStore(1,:)));
        stackMat = ccFullStore.*phaseMat;
        stackTrace = mean(stackMat,2);
        %[PWSTrace(:,ccStack),tVec] = PWS(ccFullStore,linearStack(:,ccStack),fSamp,K,nuStack,lagVal);
        ccStoreFinal(:,ccStack) = stackTrace;
        ccStack = ccStack+1;
        ccFullStore = [];
    end
end