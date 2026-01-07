% group velocity correction code, assuming the noise can be composed of a
% set of plane waves
% S. Koley, 2022

clear;
%close all;
%fPath = ['B:\LimburgBigSurvey1CC-Pair\TestSubArrays\'];
fPath = ['B:\LimburgBigSurvey1CC-Pair\SubArrayCC\'];
fName = 'subArray05FullNewCC.mat';
load([fPath,fName]);
% all the cross-correlations along with attributes are loaded now
% load the theoretical group vel from Terziet small survey
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\DataAnalysisPassive2\');
load('thGrpVelSmooth.mat');
% first get the angle of all the station pairs with respect to +ve
% x-axis, values range from 0-360 degrees
nPairs = length(rayAttribute(:,1));
for i = 1:1:nPairs
    rayAngle(i,1) = getLineAngle(rayAttribute(i,1:2),rayAttribute(i,3:4));
end
rayDist = rayAttribute(:,5);

velVec = 150:5:1000;
theta = 0:2:360;
fExt = (2.0:0.2:8)';
nF = length(fExt);
alpha = 100; beta = 3;
fSamp = 25; % downsampled in the python code
lenCC = length(ccStoreFinal(:,1));
% now get the envelopes, but for every frequency band of interest
% first do fft of the entire data set
[fftCC] = fft(ccStoreFinal,lenCC,1);
fCC = fSamp*linspace(0,1,lenCC);

% set the frequencies greater than nyquist freq, = 0
signFCC = repmat(fCC',1,nPairs);
SA = (1+sign(fCC))'.*fftCC;
SA((fCC>fSamp/2),:) = 0;
% 
tIFFT = ((-(lenCC-1)/2):1:((lenCC-1)/2))/fSamp;
tNew = tIFFT;
%tNew = tIFFT(100:900);
% %now filter the amplitude spectra and transform it back to time
for fNo = 1:1:length(fExt)
    [filterVal] = getGaussFilter(fExt(fNo,1),alpha,fCC,beta);
    %for pairNo = 1:1:nPairs
    for pairNo = 1:1:nPairs
        SAFilt(:,pairNo) = SA(:,pairNo).*filterVal';
        SAIFFT(:,pairNo) = abs((ifft(SAFilt(:,pairNo))));
        SAIFFT(:,pairNo) = SAIFFT(:,pairNo)/max(SAIFFT(:,pairNo));
    end
    %SAIFFTFull(:,fNo) = SAIFFT(100:900,6);
    %SAIFFT = SAIFFT(100:900,:);
    fCCNew = fSamp*linspace(0,1,lenCC);   
    
    beamAdd = zeros(length(velVec),length(theta));
    for i = 1:1:length(velVec)
        for j = 1:1:length(theta)
            for k = 1:1:nPairs
                tS = rayDist(k,1)*cosd(theta(j)-rayAngle(k,1))/velVec(i);
                tSInd = find(tNew>=tS,1,'first');
                beamNow = SAIFFT(tSInd,k);
                beamAdd(i,j) = beamAdd(i,j) + beamNow;
            end
        end
    end
    SAIFFT = [];
    maxValue = max(beamAdd(:));
    [r,c] = find(beamAdd==maxValue);
    velNow(fNo,1) = velVec(r(1,1)); azimuthNow(fNo,1) = theta(c(1,1));
    disp(['Max vel is ',num2str(velNow(fNo,1)),' and azimuth is ',num2str(azimuthNow(fNo,1))]);
    %plot(fExt(fNo,1),velNow);
    
    figure(10)
    surf(theta,velVec,beamAdd);
    shading interp;
    colorbar;
    view(2)
    disp('One freq done');
end
figure(11)
subplot(1,2,1)
hold on;
plot(fExt,velNow);
plot(2.6:0.25:8,vgBC,'r--'); 
plot(2.6:0.25:8,1.35*vgBC,'k--');
plot(2.6:0.25:8,0.65*vgBC,'r--');
hold off;

subplot(1,2,2)
hold on;
plot(fExt,azimuthNow);
hold off;
