% test FTAN
%load('NorthStack.mat');
close all;
load('thGrpVelSmooth.mat');
%minAlpha = 300; maxAlpha = 500;
minAlpha = 30; maxAlpha = 50;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1:0.25:3)';
plotSet = 1;
lagVal = -500:1:500;

caseString = ["causal";"acausal";"symm"];
ccPair = ccStore(:,125);
distNow = 2700;
%ccPair = ccStoreFinal;
%distNow = rayAttribute(1,5);

vMin = 1; vMax = 3; vInterval = 0.05;
%[w] = vel_taper2(tVec,1/fSamp,distNow/1000,vMin, vMax, vInterval );
%ccPair = ccPair.*(w');
%[newCC,tt] = doVelTaper(ccPair,fSamp,tArray,vMin,vMax,vInterval,distNow);
newCC = ccPair;
for signalType = 1:1:length(caseString)
    caseType = caseString(signalType,1);
    for filterNo = 1:1:length(alphaVal)
        plotNum = (signalType-1)*3+1;
        
        [fExt,distOut,tPU,velPickUpdated,phPick,spectSNR(:,signalType),figOut] = ...
                  FTANNewUpd(newCC,lagVal,minAlpha,maxAlpha,fSamp,fExt,distNow,...
                   plotSet,caseType,plotNum,Beta);
        %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
        disp(['Distance is = ', num2str(distOut)]);
        
        figure(104)
        hold on;
        plot(fExt,velPickUpdated,'-o','color',colorsAll(signalType,:));
        %plot(fExt,velPickUpdated,'color',colorsAll(signalType,:));
        plot(2.6:0.25:8,vgBC,'cyan','LineWidth',2)
        title(['Case = ', caseString(signalType,1)]);
        %figure(104);hold on;
        ylim([100,3000]);
        xlim([1,3]);
        
        figure(105);
        hold on;
        plot(fExt,spectSNR(:,signalType),'-o','color',colorsAll(signalType,:));
        plot(fExt,5*ones(length(fExt),1),'b');
        set(gca,'YScale','log');
    end
    
end
