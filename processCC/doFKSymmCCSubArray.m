% this script performs the FK transform of the symetric correlations

clear; close all;

% load all the cross-correlations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray0511.mat');
%load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCs.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray09vFKNew.mat');
kValPick = vAll(:,1)./vAll(:,2);
%highVel = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\StarsJan\OnlyFund\Star1\PhMode1FK.txt');
highVel = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\OvertFKSmooth.txt');
kValHigh = highVel(:,1)./highVel(:,2);

%load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCs.mat');
minDist = min(rayAttribute(:,5));
maxDist = max(rayAttribute(:,5));

distBin = minDist:25:maxDist;
nBin = length(distBin);
load('HdHp1Hz.mat');
load('HdHp3Hz.mat');
load('HdBp1_5to3Hz.mat');
%% plot the subArray
% load the node locations
load('nodeLocationsCartesian.mat');

% get the location of the subarray nodes
for i = 1:1:length(stnList)
    stnInd = find(allStn==stnList(i));
    stnLocUse(i,1:2) = nodeLocationsCartesian(stnInd,2:3);
end

figure(21);hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo','MarkerSize',6);
plot(stnLocUse(:,1),stnLocUse(:,2),'ro','MarkerSize',6,'MarkerFaceColor','r');
hold off;
%%
for i = 1:1:(nBin-1)
    distInd = find(rayAttribute(:,5)>distBin(i) & rayAttribute(:,5)<distBin(i+1));
    ccNow = ccStoreFinal(:,distInd);
    ccDistAvg = mean(ccNow,2);
    % make it symmetric;
    %ccDistAvg = filtfilt(HdBp1_5to3Hz.Numerator,1,ccDistAvg);
    ccMat(:,i) = [ccDistAvg(501,1);ccDistAvg(502:end,1)+flipud(ccDistAvg(1:500,1))];
    ccMat(:,i) = filtfilt(HdHp1Hz.Numerator,1,ccMat(:,i));
%     ccMat(:,i) = filtfilt(HdHp3Hz.Numerator,1,ccMat(:,i));
    
    ccMat(:,i) = ccMat(:,i)/max(abs(ccMat(:,i)));
end

% check for NaN and replace with zero
for i = 1:1:(nBin-1)
    if(isnan(ccMat(:,i)))
        ccMat(:,i) = 0;
    end
end
%%
figure(1);
xArray = (distBin(1:end-1)+distBin(2:end))/2;
plotseis(ccMat,tArray(501:end),xArray',1);

figure(3);
imagesc(xArray,tArray(501:end),ccMat);
xlabel(['Station pair offset (m)']);
ylabel(['Time lag (s)']);
ylim([0,10]);
colorbar; colormap('gray');

%%
[spec,f,kx]=fktran(ccMat,tArray(501:end),xArray);
for i = 1:1:length(spec(:,i))
    specAbs(i,:) = abs(spec(i,:))/max(abs(spec(i,:)));
end

%%
figure(2); hold on;
surf(kx,f,specAbs);
shading interp;
view(2);
colorbar;
xlim([0,0.0075]);
ylim([0.5,5]);
xlabel('wavenumber (1/m)');
ylabel('Frequency (Hz)');
hold on;
plot3(kValPick,vAll(:,1),10*ones(length(vAll),1),'ko','MarkerSize',4);
plot3(kValHigh,highVel(:,1),10*ones(length(highVel),1),'ko','MarkerSize',4);

hold off;
%%
% get the peaks from the FK transform itself
kPos = kx(kx>=0);
fSt = 1; fEnd = 5;
fStartInd = find(f>=fSt,1,'first');
fEndInd = find(f>=fEnd,1,'first');

fC = 1;
for fNo = fStartInd:fEndInd
    specPos = specAbs(fNo,kx>=0);
    [pkVal,pkLoc] = findpeaks(specPos);
    kStruct(fC).val = kPos(pkLoc);
    fC = fC+1;
end

%%
figure(4); hold on;
surf(kx,f,specAbs);
shading interp;
view(2);
colorbar;
xlim([0,0.0075]);
ylim([0.5,5]);
xlabel('wavenumber (1/m)');
ylabel('Frequency (Hz)');

fC = 1;
for fNo = fStartInd:fEndInd
    nVal = length(kStruct(fC).val);
    plot3(kStruct(fC).val,f(fNo)*ones(nVal,1),1000*ones(nVal,1),'ko','MarkerSize',6);
    fC = fC+1;
end
hold off;
%%
% also now show the nodes in the subArray with altitude
A = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
refStn = "WUJ2A";
refLong = 5.9060966; refLat = 50.7513927;

refLatLong = [refLat,refLong];
R = 6371000; % in meters
for i = 1:1:length(allStn)
    endStnInd = find(allStn==allStn(i,1),1,'first');
    latLongNow = A(endStnInd,1:2);
    [xlong,ylat] = calculatedist(latLongNow,refLatLong,R);
    nodeLocCart(i,1:3) = [xlong,ylat,A(endStnInd,3)];
end

R = 6371000; % in meters
for i = 1:1:length(stnList)
    subArrayStnInd = find(allStn==stnList(i,1),1,'first');
    latLongNow = A(subArrayStnInd,1:2);
    [xlong,ylat] = calculatedist(latLongNow,refLatLong,R);
    nodeLocCartSubArray(i,1:2) = [xlong,ylat];
end

xAlt = min(nodeLocCart(:,1)):200:max(nodeLocCart(:,1));
yAlt = min(nodeLocCart(:,2)):200:max(nodeLocCart(:,2));
[xyAltA,xyAltB] = meshgrid(xAlt,yAlt);
% xyAltA = xyAltA(:);
% xyAltB = xyAltB(:);

% now get the altitude map
altIntp = griddata(nodeLocCart(:,1),nodeLocCart(:,2),nodeLocCart(:,3),...
    xyAltA,xyAltB,'cubic');
altIntpMat = vec2mat(altIntp,length(xAlt));

figure(41)
h = imagesc(xAlt,yAlt,altIntpMat);
set(h, 'AlphaData', ~isnan(altIntpMat))
set(gca,'YDir','normal')
colorbar; colormap('jet');
hold on;
plot(nodeLocCart(:,1),nodeLocCart(:,2),'mo','MarkerSize',8,'MarkerFaceColor','k');
plot(nodeLocCartSubArray(:,1),nodeLocCartSubArray(:,2),'ko','MarkerSize',8,'MarkerFaceColor','r');
hold off;