% this script was written to get the virtual noise gather

clear; close all;

fPathBase = 'D:\LimburgSurvey1CC\stack\';

% load all the station locations, which is the metadata
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
load('nodeLocationsCartesian');
nStn = length(allStn);
refStationChar = 'XBKUA';
refStationInd = find(allStn==string(refStationChar));
refLoc = nodeLocationsCartesian(refStationInd,2:3);
refLocX = refLoc(1,1); refLocY = refLoc(1,2);

R = 6371000; % radius of earth in meters

remStn = setdiff(allStn,string(refStationChar));

A = nodeLocationsCartesian;
for i = 1:1:length(remStn)
    x = find(allStn==remStn(i,1));
    offset(i,1) = sqrt((A(x,2)-refLocX)^2 + (A(x,3)-refLocY)^2);
end

% now loop through days
% loop through hours of the day

virtMat = zeros(1001,length(remStn));
ccCount = zeros(length(remStn),1);

fPathRef = [fPathBase,refStationChar,'.mat'];
if(exist(fPathRef))
    load(fPathRef);
    % check how many stations are present
    z = string(stnStore);
    stnEndRem = setdiff(remStn,z);    
    % loop across and fill up the matrix
    for i = 1:1:length(z)
        a = find(remStn==z(i,1));
        if(~isempty(a))
            virtMat(:,a) = virtMat(:,a) + ccStoreFinal(:,i);
            ccCount(a,1) = ccCount(a,1)+1;
        end
    end
else
    % base file is absent so you have to open all remaining files
    stnEndRem = remStn;
end
% try to open remaining files

for i = 1:1:length(stnEndRem)
    fPathRef = [fPathBase,char(stnEndRem(i,1)),'.mat'];
    if(exist(fPathRef))
        % load the file
        load(fPathRef);
        z = string(stnStore);
        stnEndRemInd = find(z==string(refStationChar));
        if(~isempty(stnEndRemInd))
            fillInd = find(remStn==stnEndRem(i,1));
            virtMat(:,fillInd) = flipud(ccStoreFinal(:,stnEndRemInd));
            ccCount(fillInd,1) = ccCount(fillInd,1)+1;
        end
    else
        % thaat pair is absent for the hour of the day
    end
end

load('HdBp1_3Hz.mat');
for i = 1:1:length(remStn)
    virtMat(:,i) = filtfilt(HdBp1_3Hz.Numerator,1,virtMat(:,i));
    virtMat(:,i) = virtMat(:,i)/max(abs(virtMat(:,i)));
end


tArray = (-500:1:500)/25;
figure(1)
plotseis(virtMat,tArray,offset/1000,1);

figure(2);
imagesc(offset/1000,tArray,virtMat);
