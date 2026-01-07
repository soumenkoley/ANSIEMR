% this script was written to pick the travel times between station pairs by
% using the theoretical travel times as a reference
%clear; close all;

addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');
% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-corr and the theoretical travel times
load('grpTheoTimes.mat');
% also load the vgBC from Terziet circular array campaign
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\vgBCTerziet.mat');
% load the group velocity picks from initial analysis
grpCurve = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\Star62\FundUpd.txt');

stnA = "YEL9A"; stnB = "X4MBA";

% FTAN parameters
minAlpha = 30; maxAlpha = 50;
%minAlpha = 1000; maxAlpha = 1200;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.5:0.05:10)';
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

stnAInd = find(allStn==stnA);
stnBInd = find(allStn==stnB);

stnANumInd = find(rayAttribute(:,7)==stnAInd);
stnBNumInd = find(rayAttribute(:,8)==stnBInd);

indAB = intersect(stnANumInd,stnBNumInd,'stable');
if(isempty(indAB))
    stnANumInd = find(rayAttribute(:,7)==stnBInd);
    stnBNumInd = find(rayAttribute(:,8)==stnAInd);
    indAB = intersect(stnANumInd,stnBNumInd,'stable');
    ccNow = flipud(ccStoreFinal(:,indAB));
    stnDist = rayAttribute(indAB,5);
else
    ccNow = ccStoreFinal(:,indAB);
    stnDist = rayAttribute(indAB,5);
end

%stnDist = 900;
%ccNow = ccNowT;
grpTimeIni = stnDist./grpCurve(:,2);
grpTimeUp = stnDist./(1.3*grpCurve(:,2));
grpTimeLow = stnDist./(0.7*grpCurve(:,2));

% now proceed with FTAN of ccNow
% first get the theoretical time
thTimeNow = thTime(:,indAB);
fTh = 1.3:0.1:2.5;
% interpolate this time on fExt vector
thTimeNowIntp = interp1(fTh,thTimeNow,fExt);
%figure(316)
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\PWSStackTest\azSum1500.mat');
%ccNow = ccTest;
for s = 1:1:length(caseString)
    caseType = caseString(s,1);
    
    plotNum = (s-1)*3+1;
    
    [fExt,distOut,SAOut,tOut,symmCC,figOut] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
    %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
    
    disp(['Distance is = ', num2str(distOut)]);
    
    % now the function to identify the arrival times
    [tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);
    
    % now get the SNR corresponding to these time-picks
    
    
    spectralSNR(:,s) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPick(:,s));
    
    %subplot(1,3,s)
    figure(200+s)
    surf(fExt,tOut,SAOut);
    shading interp;
    colorbar; colormap('jet');
    title(caseString(s,1));
    ylim([0,6]);
    view(2);
    hold on;
    %plot3(fExt,tPick(:,s),1000*ones(length(fExt),1),'ko','MarkerFaceColor','k');
    %plot3(1.3:0.1:2.5,thTimeNow,1000*ones(length(thTimeNow),1),'m-','LineWidth',3);
    %plot3(grpCurve(:,1),grpTimeIni,2000*ones(length(grpCurve(:,1)),1),'g-','LineWidth',3);
    %plot3(grpCurve(:,1),grpTimeUp,2000*ones(length(grpCurve(:,1)),1),'g--','LineWidth',3);
    %plot3(grpCurve(:,1),grpTimeLow,2000*ones(length(grpCurve(:,1)),1),'g--','LineWidth',3);
    
    %plot3(vgBC(:,1),stnDist./vgBC(:,2),2000*ones(22,1),'m--');
    set(gca,'YDir','reverse');
    hold off;
    
end
% now time to use the most closest travel time
for i = 1:1:length(fExt)
    % get the minimum error of the three cases, symm, acausal and causal
    b = (~isnan(tPick(i,:)));
    nanFlag  = sum(b);
    if(nanFlag>0)
        % means there is at least one not-Nan element
        [~,minInd] = min(tErrOut(i,:));
        tFinalPick(i,1) = tPick(i,minInd);
        tFinalPick(i,1) = (tFinalPick(i,1)+thTimeNowIntp(i,1))/2;
        tFinalError(i,1) = abs(tFinalPick(i,1)-thTimeNowIntp(i,1));
        finalSNR(i,1) = spectralSNR(i,minInd);
    else
        tFinalPick(i,1) = NaN;
        tFinalError(i,1)= NaN;
        finalSNR(i,1) = NaN;
    end
end
% now smooth it further
%tFinalPick = smooth(tFinalPick);

figure(317)
subplot(2,1,1)
hold on;
plot(fExt,thTimeNowIntp,'r');
plot(fExt,tFinalPick,'b*');
hold off;

subplot(2,1,2)
plot(fExt,finalSNR,'*');

figure(318)
hold on;
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ko','MarkerSize',6);
%hold on;
raySt = [nodeLocationsCartesian(stnAInd,2),nodeLocationsCartesian(stnAInd,3)];
rayEnd = [nodeLocationsCartesian(stnBInd,2),nodeLocationsCartesian(stnBInd,3)];
plot(nodeLocationsCartesian(stnAInd,2),nodeLocationsCartesian(stnAInd,3),'ko','MarkerFaceColor','k','MarkerSize',6);
plot(nodeLocationsCartesian(stnBInd,2),nodeLocationsCartesian(stnBInd,3),'ko','MarkerFaceColor','k','MarkerSize',6);
plot([raySt(1,1),rayEnd(1,1)],[raySt(1,2),rayEnd(1,2)],'k-');

ccSym = [ccNow(501,1);ccNow(502:1001,1)+flipud(ccNow(1:500))];

% execute the new version of the code

% caseType = caseString(3,1);
% [fExt,distOut,SAOut,tOut,symmCC,figOut] = FTANPhVel(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
%         plotSet,caseType,plotNum,Beta);
%     
fftCC = fft(ccSym);
freqCC = linspace(0,1,501)*25;
figure(245)
subplot(2,1,1)
plot(freqCC,abs(fftCC));
subplot(2,1,2)
plot(freqCC,unwrap(angle(fftCC)));

fStart = 1; fEnd = 10;
fStartInd = find(freqCC>=fStart,1,'first');
fEndInd = find(freqCC>=fEnd,1,'first');
tPhase = unwrap(angle(fftCC(fStartInd:fEndInd,1)))./(2*pi*(freqCC(fStartInd:fEndInd))');
vPhase = distOut./tPhase;

%figure(316)
%subplot(1,3,3); hold on;
%plot3(freqCC(fStartInd:fEndInd),abs(tPhase),1000*ones(length(tPhase),1),'k','LineWidth',2);

figure(98)
hold on;
plot(freqCC(fStartInd:fEndInd),abs(vPhase))