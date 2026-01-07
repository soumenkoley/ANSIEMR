%testPlot
% also do ftan of noise gather
close all;
ccStoreFilt = [];
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\subArray05\';
A = load([fPath,'Overt.txt']);
load('HdBp2_4Hz.mat')
%now filter and see which correlations are symmetric
ccStoreFilt = ccStoreFinal;
for i = 1:1:length(ccStoreFinal(1,:))
    ccStoreFilt(:,i) = filtfilt(HdBp2_4Hz.Numerator,1,ccStoreFinal(:,i));
end

for i = 1:1:length(ccStoreFilt(1,:))
    ccStoreFilt(:,i) = ccStoreFilt(:,i)/max(abs(ccStoreFilt(:,i)));
end

% for i = 1:1:length(ccStoreFinal(1,:))
%     tapOut = vel_taper2(tArray,1/fSamp,rayAttribute(i,5)/1000,0.15,1.2,0.05);
%     ccStoreFilt(:,i) = ccStoreFilt(:,i).*tapOut';
% end

load('HdBp1_5Hz.mat')
%now filter and see which correlations are symmetric
for i = 1:1:length(ccStoreFinal(1,:))
    ccStoreFilt(:,i) = filtfilt(HdBp1_5Hz.Numerator,1,ccStoreFinal(:,i));
end

for i = 1:1:length(ccStoreFilt(1,:))
    ccStoreFilt(:,i) = ccStoreFilt(:,i)/max(abs(ccStoreFilt(:,i)));
end

% for i = 1:1:length(ccStoreFinal(1,:))
%     tapOut = vel_taper2(tArray,1/fSamp,rayAttribute(i,5)/1000,0.15,1.2,0.05);
%     ccStoreFilt(:,i) = ccStoreFilt(:,i).*tapOut';
% end

% now sum the traces in distance bins
dx = 50;
distBin = min(rayAttribute(:,5)):dx:max(rayAttribute(:,5));
distBin = [distBin,max(rayAttribute(:,5))];
nBins = length(distBin(1,:));
stackedTrace = [];

vMin = 0.3; vMax = 0.8; vInterval = 0.05;
fSamp = 25;
for i = 1:1:(nBins-1)
    [ind,~] = find(rayAttribute(:,5) >= distBin(i) & rayAttribute(:,5)<distBin(i+1));
    stackedTrace(:,i) = mean(ccStoreFilt(:,ind),2,'omitnan');
    stackedTrace(:,i) = stackedTrace(:,i)/max(abs(stackedTrace(:,i)));
    distNow = (distBin(i) + distBin(i+1))/2;
    
    %[w] = vel_taper2(tArray,1/fSamp,distNow/1000,vMin, vMax, vInterval );
    %stackedTrace(:,i) = stackedTrace(:,i).*(w');
end
% however there can be NaN for some distance bins, so interpolate
newStacked = stackedTrace(:,~isnan(stackedTrace(1,:)));
newDistBin = distBin(1,~isnan(stackedTrace(1,:)));

figure(30);
plotseis(newStacked,tArray,newDistBin,1);

interpDist = min(newDistBin(1,1)):dx:max(newDistBin(1,end));
interpStack = [];
for i = 1:1:length(tArray)
    interpStack(i,:) = interp1(newDistBin,newStacked(i,:),interpDist);
end

figure(31);
plotseis(interpStack,tArray,interpDist,1);
% now symmetrize it

figure(39);
imagesc(interpDist,tArray,interpStack);

[spec,f,kx]=fktran(interpStack,tArray,interpDist);
% normalize every frequency bin by maximum
fStart = 1; fEnd = 5;
fStartInd = find(f>=fStart,1,'first');
fEndInd = find(f>=fEnd,1,'first');
fSmall = f(fStartInd:fEndInd,1);
for i = fStartInd:1:fEndInd
    spec(i,:) = abs(spec(i,:))/max(abs(spec(i,:)));
end
spec = spec(fStartInd:fEndInd,:);
figure(32)
surf(kx,f(fStartInd:fEndInd),abs(spec));
shading interp;
% plot the loaded overtone on this
for i = 1:1:length(A)
    kObs(i,1) = A(i,1)/A(i,2);
end
hold on;
plot3(kObs,A(:,1),1000*ones(length(kObs),1),'k--','LineWidth',1);
plot3(-kObs,A(:,1),1000*ones(length(kObs),1),'k--','LineWidth',1);
view(2)
ylim([0,5]);
xlim([-0.005,0.005]);

% sum f-k spectrum along positive and negative kx
k0Ind = find(kx==0);
specNew=[];
for i = 1:1:length(spec(:,1))
    specNew(i,:) = [spec(i,k0Ind),(spec(i,(k0Ind+1):end)+...
        fliplr(spec(i,2:(k0Ind-1))))/2,spec(i,1)];
end

figure(35)
surf(-fliplr(kx(1:k0Ind)),fSmall,abs(specNew));
shading interp;
view(2)