% this script was writte to assemble the cross-correlations generated from
% the python code getCrossCorr.ipynb
% we generate only cross-corr pairs

clear; %close all;
dayVec1 = [13:1:30]; dayVec2 = [1:6];
doy = 318:1:341;
dayVec = [dayVec1,dayVec2];
monthVec = [11*ones(1,length(dayVec1)),12*ones(1,length(dayVec2))];
yearVec = [20*ones(1,length(dayVec1)),20*ones(1,length(dayVec2))];

hrVec = 0:1:23;

dataPath = 'D:\LimburgSurvey1CC';

stnA = 'YRNVA'; stnB = 'YHNNA';

% Also read all the stan names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');

ccAvg = zeros(1001,1);
ccFullStore = [];
fSamp = 25;
tArray = [(-500:1:-1),(0:1:500)]/fSamp;
% now loop
ccCount = 1;
for dayNo = 1:1:length(dayVec)
    if(dayVec(dayNo)<10)
        dayStr = ['0',num2str(dayVec(dayNo))];
    else
        dayStr = num2str(dayVec(dayNo));
    end
    
    %foldPathA = [num2str(yearVec(dayNo)),num2str(monthVec(dayNo)),dayStr];
    foldPathA = ['Day',num2str(doy(dayNo))];
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
            
            ccAvg = ccAvg + ccTemp;
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
% try phase weighted stacking
%K = 6; nuStack = 2;
%lagVal = (-500:1:500);

load('HdBp1_5to3Hz.mat');
for i = 1:1:length(ccFullStore(1,:))
    ccFullStore(:,i) = filtfilt(HdBp1_5to3Hz.Numerator,1,ccFullStore(:,i));
end

for i = 1:1:length(ccFullStore(1,:))
    ccFullStore(:,i) = ccFullStore(:,i)/max(abs(ccFullStore(:,i)));
end

linearStack = mean(ccFullStore,2);
%[PWSTrace,tVec] = PWS(ccFullStore,linearStack,fSamp,K,nuStack,lagVal);

