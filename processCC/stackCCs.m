% this script was written to stack the CCs over all hours, days
clear; close all;
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');
dayVec1 = 12:1:30; dayVec2 = [1:1:6];
dayNumYear = 317:1:341;
dayVec = [dayVec1,dayVec2];
monthVec = [11*ones(1,length(dayVec1)),12*ones(1,length(dayVec2))];
yearVec = [20*ones(1,length(dayVec1)),20*ones(1,length(dayVec2))];

hrVec = [0:23];

%dataPath = 'B:\LimburgBigSurvey1CC-Pair\';
dataPath = 'D:\LimburgSurvey1CC\';
storePath = 'D:\LimburgSurvey1CC\stack\';
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
stnLocs = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
refLong = 5.9060966; refLat = 50.7513927;

nStn = length(allStn);

for i = 1:1:nStn
    stnNow = allStn(i,1);
    remStn = setdiff(allStn,stnNow);
    ccStoreFinal = zeros(1001,length(remStn));
    ccCount = zeros(1,length(remStn));
    for dayNo = 1:1:length(dayNumYear)
        for hrNo = 1:1:length(hrVec)
            dayStr = ['Day',num2str(dayNumYear(dayNo))];
            hrStr = ['Hr',num2str(hrVec(hrNo))];
            stnStr = [char(allStn(i,1)),'.mat'];
            fPathFull = [dataPath,dayStr,'\',hrStr,'\',stnStr];
            
            if(exist(fPathFull))
                load(fPathFull);
                stnEnd = string(stnEnd);
                for m = 1:1:length(stnEnd)
                    stnInd  = find(remStn==stnEnd(m,1));
                    ccStoreFinal(:,stnInd) = ccStoreFinal(:,stnInd) + ccStore(:,m);
                    ccCount(1,stnInd) = ccCount(1,stnInd)+1;
                end
            end
        end
        %disp('One day')
    end
    %disp('One station!');
    ccStoreFinal = ccStoreFinal(:,(ccStoreFinal(1,:)~=0));
    stnStore = remStn((ccCount(1,:)~=0),1);
    ccCount = ccCount(:,(ccCount(1,:)~=0));
    
    
    if(~isempty(ccStoreFinal))
        % average it
        ccStoreFinal = ccStoreFinal./repmat(ccCount,1001,1);
        storeName = [char(stnNow),'.mat'];
        save([storePath,storeName],'ccStoreFinal','stnStore');
    end
    disp(['station ',num2str(i), ' done!']);
end