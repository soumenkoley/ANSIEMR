% this script was written to get virtual noise gathers with respect to one
% reference station

clear; %close all;
mainPath = 'B:\LimburgBigSurvey1CC-Pair\';

dateVec1 = [10:1:30]; dateVec2 = [1:1:2];
monthVec1 = 11*ones(1,length(dateVec1));
monthVec2 = 12*ones(1,length(dateVec2));
dateVec = [dateVec1,dateVec2];
monthVec = [monthVec1,monthVec2];
yrVec = 20*ones(1,length(dateVec));
hrVec = 0:1:23;

% create the dayStr array
for i = 1:1:length(dateVec)
    if(dateVec(i)<10)
        dayStr(i,:) = ['0',num2str(dateVec(i))];
    else
        dayStr(i,:) = num2str(dateVec(i));
    end
end

% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
refStn = "WUJ2A";
refStnInd = find(allStn==refStn,1,'first');
endStns = [allStn(1:(refStnInd-1),:);allStn((refStnInd+1):end,:)];

% define the zero matrices
ccAvg = zeros(1001,length(endStns));
ccCount = zeros(1,length(endStns));
for i = 1:1:length(dateVec)
    for j = 1:1:length(hrVec)
        subPath = [num2str(yrVec(i)),num2str(monthVec(i)),dayStr(i,:),'\Hr',num2str(hrVec(j)),'\'];
        refStnPath = [subPath,char(refStn),'.mat'];
        fullRefPath = [mainPath,refStnPath];
        if(exist(fullRefPath))
            disp(fullRefPath);
            % if reference path does not exist there is no else
            load(fullRefPath);
            %disp('loaded');
            stnEndStr = string(stnEnd);
            % now loop over endStns
            for endStnNo = 1:1:length(endStns)
                endStnInd = find(stnEndStr==endStns(endStnNo,1));
                if(~isempty(endStnInd))
                    ccAvg(:,endStnNo) = ccAvg(:,endStnNo) + ccStore(:,endStnInd);
                    ccCount(1,endStnNo) = ccCount(1,endStnNo) + 1;
                else
                    % read another file with the endstnName
                    endStnNow = char(endStns(endStnNo,1));
                    
                    [ccOut,ccFlagOut] = readEndStn(mainPath,subPath,endStnNow,refStn);
                    ccAvg(:,endStnNo) = ccAvg(:,endStnNo) + ccOut;
                    ccCount(1,endStnNo) = ccCount(1,endStnNo) + ccFlagOut;
                end
                
            end
        end
        %disp('One hour read');
    end
    disp('All hours read in day')
end
disp('Averaging');
for i = 1:1:length(ccAvg(1,:))
    if(ccCount(1,i)~=0)
        ccAvg(:,i) = ccAvg(:,i)/ccCount(1,i);
    end
end


