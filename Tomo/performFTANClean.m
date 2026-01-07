% this script was written to test the cross-correlation cleaning that I
% thought in frequency-time domain
clear; close all;

%% addpath to some folders for FTAN and other routines
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');

% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');

% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

% load the rayAttribute, cross-correlations and station names
% contain all station pairs and the start and end station coordinates
% for that pair
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCFundUpdApril.mat');

% this is a guiding velocity map for the process to work
% we pick group velocities for several station pairs
% these group velocity picks are quite reliable, and the midpoint along
% that ray is assigned that group velocity
% then Gaussian smoothing is applied by using the group velocity at those
% points and approximate group velocity maps are created
% initial group travel times between all station pairs are generated on
% this map by using straight ray travel times
% these travel times between all pairs are treated as guiding light
fTh = 1.0:0.05:5;
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTimeFundApril.mat');
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelOvertLimits.mat');

% load the velocity search band
%load('vBand.mat');
% this exports two variables vLow and vUp each of size nFreq x 2
% first column is frequency and second column is the velocity
load('vBandJuly.mat');
% loads an approximate band for picking group travel times
% from FTAN of raw cross-correlations
load('vBandFundAug.mat');

% the name of two stations in question
stnA = "0RPXA"; stnB = "Z3O2A";
NPi = 0;

% load the star velocity, these stars are basically the midpoints of the
% rays that was used to generate the Gaussian smoothed group velocity maps
% which serve as an initial model to work with
grpStar = load('A:\TestInver\StarsApril\Star19\FundUpdApril.txt');
% loads a variable vAll, 20x2, 1.6-3.5 Hz, phase velocity picks from FK
% analysis of a subArray of cross-correlations
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray06vFK.mat');
% loads phase velocity which were obatined from cross-correlation
% beamforming of croos-correlation pairs of a subArray
vPhCCBF = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray06vCCBF.txt');
% loads the group velocity derived from the phase velocity from FK analysis
% of cross-correlations in a subArray
vGrpArray = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray06vGrpFK.txt');

%% FTAN parametersP% FTAN parameters
minAlpha = 30; maxAlpha = 50;
Beta = 3;
colorsAll = distinguishable_colors(3);
% sampling frequency of cross-correlations
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;

% frequencies at which FTAN will be performed for the particular
% cross-correlation pair
fExt = (0.1:0.05:5)';
%fExt = (logspace(log10(0.5),log10(5),20))';
% this frequency vector is used for cleaning the cross-correlation
% pairs, and this technique I thought of, but we can improve it
fExtClean = (logspace(log10(0.5),log10(5),20))';
% a plotting flag
plotSet = 0;
% lagval of the cross-correlations, so time is from -20 - 20 ssince the
% sampling frequency is 25 Hz
lagVal = -500:1:500;
% three cases of what part of the correlation to analyze
caseString = ["causal";"acausal";"symm"];
% this is the frequency dependent velocity band which will be used fro
% cleaning the cross-correlations
% interpolate it on fExtClean
vBandUL(:,1) = interp1(vLow(:,1),vLow(:,2),fExtClean);
vBandUL(:,2) = interp1(vUp(:,1),vUp(:,2),fExtClean);

%% locating the station pair from among all the correlations
% find the station indices and hence the correct correlation
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
    thTimeNow = thTime(:,indAB);
else
    ccNow = ccStoreFinal(:,indAB);
    stnDist = rayAttribute(indAB,5);
    thTimeNow = thTime(:,indAB);
end

%% perform FTAN and time picking

% peform the analysis for the symmetric correlation, sume causal and
% acausal
caseType = caseString(3,1);

plotNum = 1;
 
% perform FTAN of the raw-correlations, before cleaning it of spurious
% energy
[fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
symmCCCopy = symmCC;
disp(['Distance is = ', num2str(distOut)]);

% a time picking routine by using the velocity band vBand
% note that this vBand is different from the velocity band used in cleaning
% the cross-correlations
[vPickSymm,tPickSymmRaw] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
%[vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
% now get the SNR corresponding to these time-picks
spectralSNRRaw(:,1) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymmRaw(:,1));

figure(1)
hold on;
surf(fExt,tOut,SAOut); shading interp;
view(2);colormap('jet');
plot3(fExt,tPickSymmRaw,ones(length(fExt),1),'k','LineWidth',2);
plot3(fTh,thTimeNow,ones(length(fTh),1),'m--','LineWidth',2);
%plot3(vBand(:,1),distOut./vBand(:,2),100*ones(length(vBand(:,1)),1),'m','LineWidth',2);
%plot3(vBand(:,1),distOut./vBand(:,3),100*ones(length(vBand(:,1)),1),'m','LineWidth',2,'HandleVisibility','off');
plot3(fExtClean(:,1),distOut./vBandUL(:,1),100*ones(length(vBandUL(:,1)),1),'m','LineWidth',2);
plot3(fExtClean(:,1),distOut./vBandUL(:,2),100*ones(length(vBandUL(:,1)),1),'m','LineWidth',2,'HandleVisibility','off');
ylim([0,10]);
plot3(vGrpArray(:,1),distOut./vGrpArray(:,2),1000*ones(length(vGrpArray),1),'c','LineWidth',2);
xlabel('Frequency (Hz)'); ylabel('Lag time (s)')
title('FTAN of raw cross-correlations');
legend({'FTAN image','timePick Raw','Time band for cleaning','grp time from FK analysis'})
hold off;

%% use the raw phase of the symmetric cross-correlation
% get phase velocity, this is quite different from Bensen et al, 2007
fftCC = fft(symmCC);
lSymmCC = 501;
fCC = linspace(0,1,lSymmCC)*fSamp;
fAInd = find(fCC>=0.5,1,'first');

angFFT = unwrap(angle(fftCC(fAInd:end,1)))+2*NPi*pi;
vv = 2*pi*distOut*fCC(fAInd:end)'./(angFFT);

figure(1008);
hold on;
plot(fCC(fAInd:end),abs(vv),'b')
xlim([1,5]); ylim([200,4000]);

figure(1009);
plot(fCC,(angle(fftCC)),'b');

figure(1010);
hold on;
plot(fCC,unwrap(angle(fftCC)),'b');
hold off;
ylabel('phase unwrapped');
%%
tIFFT = (0:1:(lSymmCC-1))/fSamp;
%now filter the amplitude spectra and transform it back to time
nF = length(fExtClean);

for fNo = 1:1:nF
    alpha(fNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(nF-1)*(fNo-1))-1;
end

alpha = fliplr(alpha);

%filterSum = zeros(1,lSymmCC);
figure(2);hold on;
ccClean = zeros(lSymmCC,1);
for fNo = 1:1:length(fExtClean)
    [filterVal] = getGaussFilterFull(fExtClean(fNo,1),alpha(fNo),fCC,Beta);
    %plot(fCC,filterVal,'b');
    filterValFull = [filterVal(1:251),fliplr(filterVal(2:251))];
    filtFFT = fftCC.*filterValFull';
    ifftCC = ifft(filtFFT);
    [tW] = genWW(vBandUL(fNo,1),vBandUL(fNo,2),distOut,fSamp,tIFFT);
    ccClean = ccClean + ifftCC.*(tW)';
end

plot(tIFFT,ccClean);
hold off;

ccClean = [zeros(500,1);ccClean];
[fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccClean,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
[vPickSymmClean,tPickSymmClean] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
%[vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
% now get the SNR corresponding to these time-picks
spectralSNRClean(:,1) = getSpectSNR(symmCCCopy,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymmClean(:,1));

[f, c, mask] = phaseVel(ccClean,fSamp,distOut,0.5,4,5,3);
opts = struct();
opts.vmin = vBand(:,3);
opts.vmax = vBand(:,2);
opts.smoothness_lambda = 1e7;
opts.post_sg_win = 11;
opts.post_sg_poly = 2;
opts.debug = true;
out = pick_group_curve_dpUpd((SAOut(:,19:end))', fExt(19:end,1), tIFFT, distOut, opts);
spectralSNRCleanN(:,1) = getSpectSNR(symmCCCopy,fExt,minAlpha,maxAlpha,Beta,fSamp,[tPickSymmClean(1:18,1);out.times_pick']);
%% phase speed from Bensen et al
[fExt,distOut,SAOut,SAPh,tOut,symmCC,figOutSymm] = ...
        FTANNewUsePh(ccClean,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);
% unwrap the phase
for i = 1:1:length(fExt)
    SAPhUnwrap(:,i) = (SAPh(:,i));
    
end
% get the phase at tPickSymmClean
for i = 1:1:length(fExt)
    tPickInd = find(tIFFT>=tPickSymmClean(i,1),1,'first');
    phNow(i,1) = SAPhUnwrap(tPickInd,i);
end
% now calculate the ph speed

for i = 1:1:length(fExt)
    N = -10;
    for j = 1:11
        sPh = tPickSymmClean(i,1)./distOut + ((phNow(i,1)) -pi/4+2*pi*N)./(2*pi*fExt(i,1)*distOut);
        vPh(j,i) = 1/sPh;
        N = N+1;
    end
%     figure(90);
%     hold on;
%     plot(fExt(i,1),vPh(:,i),'k*');
%     hold off;
end

%%
figure(3)
hold on;
surf(fExt,tOut,SAOut); shading interp;
view(2);colormap('jet');
plot3(fExt,tPickSymmClean,ones(length(fExt),1),'k','LineWidth',2);
plot3(fTh,thTimeNow,ones(length(fTh),1),'m--','LineWidth',2);
plot3(fTh,1/0.7*thTimeNow,ones(length(fTh),1),'g--','LineWidth',2);
plot3(fTh,1/1.3*thTimeNow,ones(length(fTh),1),'g--','LineWidth',2);
plot3(fExt(19:end),out.times_pick,ones(length(fExt(19:end)),1),'r','LineWidth',2);
ylim([0,10]);
hold off;

%%
fftCC = fft(symmCC);
angFFT = unwrap(angle(fftCC(fAInd:end,1)))+2*NPi*pi-pi/4;
vv = 2*pi*distOut*fCC(fAInd:end)'./(angFFT);
figure(1008);
hold on;
plot(fCC(fAInd:end),abs(vv),'r')
plot(vAll(:,1),vAll(:,2),'m');
plot(vPhCCBF(:,1),vPhCCBF(:,2),'g')
xlim([1.0,5]); ylim([200,4000]);
hold off;

figure(1009);
hold on;
plot(fCC,(angle(fftCC)),'r');
hold off;
ylabel('phase wrapped');
xlim([0,5]);

figure(1010);
hold on;
plot(fCC,unwrap(angle(fftCC)),'r');
hold off;
ylabel('phase unwrapped');
xlim([0,5]);
%%
% use the theoretical group time to get the theoretical phase velocity
c0 = 2500;
fSt = 1.0;
k0 = 2*pi*fSt/c0;
thGrpVel = distOut./thTimeNow;
thGrpVel = [fTh',thGrpVel];
fStInd = find(thGrpVel(:,1)==fSt,1,'first');
fAll = thGrpVel(fStInd:end,1);
k = 2*pi*cumtrapz(fAll(1:end),1./(1*thGrpVel(fStInd:end,2)));
kC = k0-k(1,1);
k = kC+k;
phVelOut = 2*pi*fAll(1:end)./k;
figure(1008); hold on;
plot(fAll,phVelOut,'k','LineWidth',2);
plot(fAll,1.25*phVelOut,'k--','LineWidth',2,'HandleVisibility','off');
plot(fAll,0.75*phVelOut,'k--','LineWidth',2,'HandleVisibility','off');
%plot(fExt,vPh,'m');
hold off;
legend({'phasePick raw','phasePick clean','FK Array','ph from CCBF','ph from thGrp'})
%% plot the bessel function

J0 = besselj(0,(2*pi*fAll(1:end,1)*distOut)./(phVelOut(:,1)));

figure(1011);
plot(fCC,real(fft(ccClean(501:1001,1))),'b');
hold on;
plot(fAll(:,1),J0,'r');
hold off;
xlim([1,5]);
%% plot the star vel, thGrpvel, and the current pick
figure(1012);
plot(thGrpVel(:,1),thGrpVel(:,2),'b','LineWidth',2);
hold on;
plot(grpStar(:,1),grpStar(:,2),'r','LineWidth',2);
plot(fExt,distOut./tPickSymmRaw,'k','LineWidth',2);
plot(fExt,distOut./tPickSymmClean,'m','LineWidth',2);
plot(vGrpArray(:,1),vGrpArray(:,2),'c','LineWidth',2);
plot(vBand(:,1),vBand(:,2),'m--','LineWidth',2);
plot(vBand(:,1),vBand(:,3),'m--','LineWidth',2);
plot(fExt(19:end),distOut./out.times_pick,'g','LineWidth',2);

legend({'th grp vel','star grp vel','tPIckSymmRaw','tPickSymmClean','vGrpFKArray'});
ylim([0,4000]);

spectralSNRAll = [fExt,spectralSNRRaw,spectralSNRClean,spectralSNRCleanN];

%% perform another stage of ftanclean
% figure(3); hold on;
% plot3(fExtClean,distOut/4000*ones(length(fExtClean),1),100*ones(length(fExtClean),1),'k','LineWidth',2);
% plot3(vBand(:,1),distOut./vBand(:,2),100*ones(length(vBand(:,1)),1),'m','LineWidth',2);
% plot3(vBand(:,1),distOut./vBand(:,3),100*ones(length(vBand(:,1)),1),'m','LineWidth',2);
% 
% ylim([0,max(10,max(thTimeNow)+1)]);
% [x,y] = getpts();
% yIntp1 = interp1(x,y,fExtClean);
% [x,y] = getpts();
% yIntp2 = interp1(x,y,fExtClean);
% 
% %vBandUL2 = [vBandUL(:,1),distOut./yIntp];
% vBandUL2 = [distOut./yIntp2,distOut./yIntp1];
% fftCCClean = fft(ccClean(501:1001,1));
% ccClean2 = zeros(lSymmCC,1);
% for fNo = 1:1:length(fExtClean)
%     [filterVal] = getGaussFilterFull(fExtClean(fNo,1),alpha(fNo),fCC,Beta);
%     %plot(fCC,filterVal,'b');
%     filterValFull = [filterVal(1:251),fliplr(filterVal(2:251))];
%     filtFFT = fftCCClean.*filterValFull';
%     ifftCC = ifft(filtFFT);
%     [tW] = genWW(vBandUL2(fNo,1),vBandUL2(fNo,2),distOut,fSamp,tIFFT);
%     ccClean2 = ccClean2 + ifftCC.*(tW)';
%     %filterSum = filterSum+filterValFull;
%     %disp('Stop');
% end
% hold off;
% 
% ccClean2 = [zeros(500,1);ccClean2];
% [fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
%         FTANNewUse(ccClean2,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
%         plotSet,caseType,plotNum,Beta);
% [vPickSymmClean2,tPickSymmClean2] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
% %[vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
% % now get the SNR corresponding to these time-picks
% spectralSNR(:,1) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymmClean2(:,1));
% 
% figure(4)
% hold on;
% surf(fExt,tOut,SAOut); shading interp;
% view(2);colormap('jet');
% plot3(fExt,tPickSymmClean2,ones(length(fExt),1),'k','LineWidth',2);
% plot3(fTh,thTimeNow,ones(length(fTh),1),'m--','LineWidth',2);
% ylim([0,max(10,max(thTimeNow)+1)]);
% hold off;
% 
% fftCC = fft(symmCC);
% angFFT = unwrap(angle(fftCC(fAInd:end,1)));
% vv = 2*pi*distOut*fCC(fAInd:end)'./(angFFT);
% figure(1008);
% hold on;
% plot(fCC(fAInd:end),abs(vv),'color',[1,1,0],'LineWidth',2);
% xlim([1.0,5]); ylim([200,4000]);
% hold off;
% 
% figure(1012);
% hold on;
% plot(fExt,distOut./tPickSymmClean2,'g','LineWidth',2);
% hold off;