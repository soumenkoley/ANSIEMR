% this script was written to do the Frequency-time-image analysis following
% Yao et al 2005 in order to extract the phase velocity curves
% doFTAN
%clear;
%% run this code after assembleFiniteCCPairs.m
close all;
% this function performs FTAN on a desired station pair
%addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive1\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive2\');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');
% load the picked phasse velocity for the subArray
fPathPhase = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\subArray11\';
AA = load([fPathPhase,'Overt.txt']);
% this file is from the small passive survey at Terziet
load('thGrpVelSmooth.mat');
load('HdBp3_8Hz.mat'); % a 3-8 Hz bandpass filter

%% parameters for FTAN
minAlpha = 100; maxAlpha = 200;
Beta = 3;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
n1 = "YYMLA"; n2= "Z2NAA";
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1:0.02:8)';
plotSet = 1;
%% main code starts
% now, from the stnList, find the index of the stations
stnAInd = find(stnList==n1); stnBInd = find(stnList==n2);
n1Ind = find((rayAttribute(:,7)==stnAInd) | (rayAttribute(:,7)==stnBInd));
n2Ind = find((rayAttribute(:,8)==stnAInd) | (rayAttribute(:,8)==stnBInd));

pairInd = intersect(n1Ind,n2Ind);

if(~isempty(pairInd))
    nA = rayAttribute(pairInd,7);
    nB = rayAttribute(pairInd,8);
    
    if(nA == stnAInd)
        ccPair = ccStoreFinal(:,pairInd);
    else
        ccPair = flipud(ccStoreFinal(:,pairInd));
    end
else
    error('Station pair not found!');
end
%%
% station pair read
figure(90)
ccSee = ccStoreFinal(:,pairInd);
%ccSee = filtfilt(HdBp3_8Hz.Numerator,1,ccStoreFinal(:,pairInd));
plot(tArray,ccSee);
xlim([-10,10]);
lagVal = -500:1:500; % + 20 s and -20 s on ach side
tVec = tArray;
distNow = (rayAttribute(pairInd,5));
disp(['Dist = ',num2str(distNow),' m']);

vMin = 0.1; vMax = 2; vInterval = 0.05;
%[w] = vel_taper2(tVec,1/fSamp,distNow/1000,vMin, vMax, vInterval );
%ccPair = ccPair.*(w');
%load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\AzAvgCC\cc110-38.mat');
%ccPair = meanCC;
slopeNow = rayAttribute(pairInd,6);

causalCC = ccPair(tVec>=0,1);
acausalCC = ccPair(tVec<0,1);

caseType = "symm";

switch caseType
    case 'symm'
        symmCC = (causalCC(2:end,1)+ flipud(acausalCC))/2;
        symmCC = [causalCC(1,1);symmCC];
    case 'acausal'
        symmCC = [causalCC(1,1);flipud(acausalCC)];
    case 'causal'
        symmCC = causalCC;
end
lenCC = length(symmCC);

fftCC = fft(symmCC);
fVec = [fSamp*linspace(0,1,lenCC)]';
tVec = [(1:1:lenCC)/fSamp]';
%% design the filter
nF = length(fExt);
for fNo = 1:1:nF
    alpha(fNo)=minAlpha+exp((log(maxAlpha-minAlpha+1))/(nF-1)*(fNo-1))-1;
end
alpha = fliplr(alpha);
outReal = [];

% define the cAxis
cAxUse = [200:10:4000]';
outRealInterp = [];
for fNo = 1:1:nF
    [filterVal] = getGaussFilter(fExt(fNo,1),alpha(fNo),fVec,Beta);
    fftCCNew = fftCC.*filterVal;
    %figure(21)
    %hold on;
    %plot(filterVal);
    % now make it symmetric
    [outReal] = makeFFTSymm(fftCCNew,fVec,fSamp);
    %outReal = outReal/max(abs(outReal));
    
    % convert the time axis to c-axis
    cAx = distNow./(tVec-1/(8*fExt(fNo,1)));
    cAx = flipud(cAx);
    outReal = flipud(outReal);

    outRealInterp(:,fNo) = interp1(cAx,outReal,cAxUse,'spline');
    outRealInterp(:,fNo) = outRealInterp(:,fNo)/max(abs(outRealInterp(:,fNo)));
    % also find the maximum
    [maxVal,maxInd] = max(outRealInterp(:,fNo));
    cout(fNo,1) = cAxUse(maxInd,1);
end

% now get the peak closest to the reference curve
newFreqStInd = find(fExt>=min(AA(:,1)),1,'first');
newFreqEndInd = find(fExt>=max(AA(:,1)),1,'first');
newFreq = fExt(newFreqStInd:newFreqEndInd,1);
% interpolate reference curve on this frequency vector
newRefDisp = interp1(AA(:,1),AA(:,2),newFreq);


% now ready to pick the max
k =1;
j =1;
velPick = [];
for i = newFreqStInd:1:newFreqEndInd
    [pks,locs] = findpeaks(outRealInterp(:,i));
    velPeaks = cAxUse(locs,1);
    % find the one that is closest to the ref dispersion
    [minVal,minInd] = min(abs(velPeaks-newRefDisp(j,1)));
    closestVel = velPeaks(minInd,1);
    % now we search for 10% deviation around reference curve
    velDev = 0.1*newRefDisp(j,1);
    if(minVal<velDev)
        % then store
        velPick(k,1) = fExt(i,1);
        velPick(k,2) = closestVel;
        k = k+1;
    end
    j = j+1;
end

figure(22)
surf(fExt,cAxUse,outRealInterp);
shading interp;
view(2);
hold on;plot3(AA(:,1),AA(:,2),1000*ones(length(AA),1),'k--');
plot3(velPick(:,1),velPick(:,2),1000*ones(length(velPick),1),'r*');
