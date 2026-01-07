% this script performs the FK transform of the symetric correlations

clear; close all;

% load all the cross-correlations
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCs.mat');

distBin = 0:25:6500;
nBin = length(distBin);
load('HdHp1Hz.mat');

for i = 1:1:(nBin-1)
    distInd = find(rayAttribute(:,5)>distBin(i) & rayAttribute(:,5)<distBin(i+1));
    ccNow = ccStoreFinal(:,distInd);
    ccDistAvg = mean(ccNow,2);
    % make it symmetric;
    ccMat(:,i) = [ccDistAvg(501,1);ccDistAvg(502:end,1)+flipud(ccDistAvg(1:500,1))];
    ccMat(:,i) = filtfilt(HdHp1Hz.Numerator,1,ccMat(:,i));
    ccMat(:,i) = ccMat(:,i)/max(abs(ccMat(:,i)));
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
% load the fundamental and so called overt velocities from beamforming
fundBeam = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\Fund.txt');
overtBeam = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\Overt.txt');

kFund = fundBeam(:,1)./fundBeam(:,2);
kOvert = overtBeam(:,1)./overtBeam(:,2);

plot3(kFund(:,1),fundBeam(:,1),1000*ones(length(fundBeam(:,1)),1),'mo','MarkerSize',6);
plot3(kOvert(:,1),overtBeam(:,1),1000*ones(length(overtBeam(:,1)),1),'mo','MarkerSize',6);
hold off;

%%
% get the peaks from the FK transform itself
kPos = kx(kx>=0);
fSt = 1; fEnd = 4;
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
    plot3(kStruct(fC).val,f(fNo)*ones(nVal,1),1000*ones(nVal,1),'mo','MarkerSize',6);
    fC = fC+1;
end
plot3(kFund(:,1),fundBeam(:,1),1000*ones(length(fundBeam(:,1)),1),'k*','MarkerSize',6);
plot3(kOvert(:,1),overtBeam(:,1),1000*ones(length(overtBeam(:,1)),1),'k*','MarkerSize',6);
hold off;

%%
figure(5); hold on;
surf(kx,f,specAbs);
shading interp;
view(2);
colorbar;
xlim([0,0.0075]);
ylim([0.5,5]);
xlabel('wavenumber (1/m)');
ylabel('Frequency (Hz)');

% load the FK picked velocities
fundFK = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\FundFKSmooth.txt');
overtFK = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\OvertFKSmooth.txt');

plot3(fundFK(:,1)./fundFK(:,2),fundFK(:,1),1000*ones(length(fundFK(:,1)),1),'mo','MarkerSize',6);
plot3(overtFK(:,1)./overtFK(:,2),overtFK(:,1),1000*ones(length(overtFK(:,1)),1),'ko','MarkerSize',6);
hold off;

%%
figure(6);
hold on;
plot(fundFK(:,1),fundFK(:,2),'r','LineWidth',2);
plot(fundBeam(:,1),fundBeam(:,2),'go','LineWidth',2);

plot(overtFK(:,1),overtFK(:,2),'b','LineWidth',2);
plot(overtBeam(:,1),overtBeam(:,2),'mo','LineWidth',2);
hold off;
ylabel('phase velocity (m/s)');
xlabel('Frequency (Hz)');
title('FK and Beamform phase velocities')
grid on;
box on;

%%
