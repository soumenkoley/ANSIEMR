% this script was written to pick the phase travel time between
% a pair of stations, the algorithm is still not clear in my head
% but I will try several things
clear; close all;

%% load input files
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr');
addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN\');
% load all the station names
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
% load the station coordinates
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
% load the rayAttribute, cross-correlations and station names
load('B:\LimburgBigSurvey1CC-Pair\SubArrayCC\allCCs.mat');
% also load the vgBC from Terziet circular array campaign
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\vgBCTerziet.mat');
% load the theoretical travel times
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpTime.mat');
%load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelOvertLimits.mat');

% load the velocity search band
load('vBand.mat');

stnA = "W3JWA"; stnB = "XDLDA";
% load the star velocity
%phStar = load('A:\TestInver\StarsApril\Star28\PhUpdAprilTheo.txt');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray03vFK.mat');

%% FTAN parametersP% FTAN parameters
minAlpha = 30; maxAlpha = 50;
%minAlpha = 1000; maxAlpha = 1200;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (0.1:0.05:5)';
fShow = 1.6:0.05:5;
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];

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

% updating this on Feb 05, 2025
%thTimeNow = stnDist./thGrpVel;
%% travel time based on the Terziet group velocity measurement

grpTimeIni(:,1) = vgBC(:,1); % frequecy Hz
grpTimeIni(:,2) = stnDist./vgBC(:,2); % travel time in secs

vDiff = 0.2; % allowing for a 20% devviation of the causal or acausal travel
% time to the symmetric travel times

%% perform FTAN and time picking

% peform the analysis for the symmetric correlation
caseType = caseString(3,1);

plotNum = 1;
    
[fExt,distOut,SAOut,tOut,symmCC,figOutSymm] = ...
        FTANNewUse(ccNow,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
        plotSet,caseType,plotNum,Beta);

disp(['Distance is = ', num2str(distOut)]);

%% updated on July 14, 2025
fftCC = fft(symmCC);
fVec = linspace(0,1,501)*25;
J0 = besselj(0,(2*pi*vAll(:,1)*distOut)./(vAll(:,2)));
figure(1007);
hold on;
plot(fVec,real(fftCC),'b');
plot(vAll(:,1),J0,'r');
hold off;
xlim([1,5]);

% now the function to identify the arrival times
%[tPick(:,s),tErrOut(:,s)] = pickTT(SAOut,fExt,tOut,thTimeNow,fTh);

[vPickSymm,tPickSymm] = pickTTSymmJan(SAOut,fExt,tOut,distOut,vBand,fSamp);
%[vPickSymm,tPickSymm] = pickTTSymmJanAdv(SAOut,fExt,fShow,tOut,distOut,thTimeNow,fSamp);
% now get the SNR corresponding to these time-picks
spectralSNR(:,1) = getSpectSNR(symmCC,fExt,minAlpha,maxAlpha,Beta,fSamp,tPickSymm(:,1));

figure(1)
hold on;
surf(fExt,tOut,SAOut); shading interp;
view(2);colormap('jet');
plot3(fExt,tPickSymm,ones(length(fExt),1),'k','LineWidth',2);
ylim([0,10]);
hold off;


vAll = [vPickSymm];
for i = 1:1:length(fExt)
    nanInd = isnan(sum(vAll(i,1)));
    if(~nanInd)
        vFinal(i,1) = fExt(i,1);
        vFinal(i,2) = vAll(i,1); % fill in with the symmetric velocity
    else
        vFinal(i,1) = fExt(i,1);
        vFinal(i,2) = NaN; % fill in with the symmetric velocity
    end
end

% now check if the these SNRs are less than 5

for i = 1:1:length(fExt)
    minSNR = min(spectralSNR(i,1));
    if(minSNR<5)
        vFinal(i,2) = NaN;
    end
end

figure(2);
hold on;

plot(fExt,vAll(:,1),'g','LineWidth',2);
plot(fExt,vAll(:,1),'go');
plot(fExt,vFinal(:,2),'k*');
set(gca,'YScale','log');
xlim([1,5]);
hold off;

% now the best fitting bessel function part
f1 = 0.5; f2 = 5;
f1Ind = find(fVec>=f1,1,'first');
f2Ind = find(fVec>=f2,1,'first');
fSmall = (fVec(f1Ind:f2Ind))';

fftCCReal = real(fftCC(f1Ind:f2Ind,1));
fftCCReal = fftCCReal./max(abs(fftCCReal));
fftCCReal = smooth(fftCCReal);

% % Initial guess: constant phase velocity (e.g., 3000 m/s everywhere)
% c0 = 1000*ones(length(fSmall),1);
% % load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray03vFK.mat');
% % 
% % c0 = interp1(vAll(:,1),vAll(:,2),fSmall,'pchip','extrap');
% 
% % Regularization parameter (tune this):
% lambda = 1; % Try 1, 10, 100 and see
% 
% % Use fminunc or other solvers:
% lb = 100 * ones(length(fSmall),1);
% ub = 3000 * ones(length(fSmall),1);
% 
% obj_fun = @(c_vec) bessel_dispersion_residual(c_vec,fSmall,distOut,fftCCReal,lambda);
% 
% options = optimoptions('lsqnonlin','Display','iter');
% 
% cOut = lsqnonlin(obj_fun, c0, lb, ub, options);

% phi = acos(fftCCReal); % Safe because Re_Corr ? [-1, 1]
% 
% %Unwrap the phase to handle multiple cycles:
% 
% n = zeros(size(fSmall));
% for i = 2:length(fSmall)
%     if fftCCReal(i-1,1)*fftCCReal(i,1)<0
%         n(i:end) = n(i:end) + 1; % Each zero crossing adds ? to phase
%     end
% end
% 
% % Reconstruct total phase
% theta = phi + n*pi;
% 
% % Invert for c(f)
% cIni = (2*pi.*fSmall.*distOut)./(theta + pi/4);
% cIni = smooth(fSmall,cIni, 0.2, 'loess'); % 20% smoothing
vv = abs(2*pi*distOut*fVec'./(unwrap(angle(fftCC))));
cIni = vv(f1Ind:f2Ind,1);

JIni = besselj(0, (2*pi.*fSmall.*distOut)./cIni);

figure(2024);
plot(fSmall,fftCCReal, 'b');
hold on;
plot(fSmall,JIni, 'r');
xlabel('Frequency (Hz)');
ylabel('Bessel value');
title('Bessel fitting');
grid on;

figure(2025);
plot(fSmall,cIni, 'b-o');
xlabel('Frequency (Hz)');
ylabel('Phase Velocity (m/s)');
title('Estimated Phase Velocity Dispersion');
grid on;


% % Regularization parameter (tune this):
lambda = 0.1; % Try 1, 10, 100 and see
% Use fminunc or other solvers:
lb = 100 * ones(length(fSmall),1);
ub = 4000 * ones(length(fSmall),1);

obj_fun = @(c_vec) bessel_dispersion_residual(c_vec,fSmall,distOut,fftCCReal,lambda);

options = optimoptions('lsqnonlin','Display','iter');

cOut = lsqnonlin(obj_fun, cIni, lb, ub, options);

JOut = besselj(0, (2*pi.*fSmall.*distOut)./cOut);

figure(2026);
plot(fSmall,fftCCReal, 'b');
hold on;
plot(fSmall,JOut, 'r');
xlabel('Frequency (Hz)');
ylabel('Bessel value');
title('Bessel fitting');
grid on;

figure(2027);
plot(fSmall,cOut, 'b-o');
xlabel('Frequency (Hz)');
ylabel('Phase Velocity (m/s)');
title('Estimated Phase Velocity Dispersion');
grid on;

function L = bessel_dispersion_loss(c_vec, f, d, Re_Corr, lambda)
    J = besselj(0, (2 * pi .* f .* d) ./ c_vec);
    misfit = norm(J - Re_Corr);
    
    % Smoothness penalty (second derivative)
    dc = diff(c_vec);
    d2c = diff(dc);
    smoothness = norm(d2c);
    
    L = misfit + lambda * smoothness;
end

function residual = bessel_dispersion_residual(c_vec, f, d, Re_Corr, lambda)
    model = besselj(0, (2 * pi .* f .* d) ./ c_vec);
    data_residual = model - Re_Corr;
    
    % Smoothness residual (second derivative)
    dc = diff(c_vec);
    d2c = diff(dc);
    smoothness_residual = sqrt(lambda) * d2c; % Scale by sqrt(lambda)
    
    % Combine into one residual vector
    residual = [data_residual; smoothness_residual];
end