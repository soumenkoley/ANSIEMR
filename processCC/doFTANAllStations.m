% this script was written to do FTAN on all station pairs
clear; close all;

addpath('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\FTAN');

%dataPath = 'B:\LimburgBigSurvey1CC-Pair\';
dataPath = 'D:\LimburgSurvey1CC\stack\';

% read the station list
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');
A = load('nodeLocationsCartesian.mat');
A = A.nodeLocationsCartesian;
nNodes = length(A(:,1));

%% FTAN parameters
minAlpha = 300; maxAlpha = 500;
Beta = 3;
alphaVal = 300;
%colorsAll = distinguishable_colors(length(alphaVal));
colorsAll = distinguishable_colors(3);
fSamp = 25; % remember you downsampled it from 250-25 Hz
dt = 1/fSamp;
fExt = (1:0.05:3)';
plotSet = 0;
lagVal = -500:1:500;
caseString = ["causal";"acausal";"symm"];
%%
% main loop starts here
ccCount = 1;
G = "0GNDA";
GInd = find(allStn==G);

for i = GInd:1:GInd
    % load the file first
    fPath = [dataPath,char(allStn(i,1)),'.mat'];
    
    if(exist(fPath))
        load(fPath);
        disp('File loaded!');
        refStnInd = find(allStn==allStn(i,1));
        refLoc = A(refStnInd,:);
        % now loop across end stations
        for j = 1:1:length(stnStore)
            eStnInd = find(allStn==stnStore(j,1));
            endLoc = A(eStnInd,:);
            rayAttNow(1:2,1) = [refLoc(1,1);endLoc(1,1)];
            rayStart = refLoc(1,2:3); rayEnd = endLoc(1,2:3);
            rayAttNow(3:4,1) = rayStart';
            rayAttNow(5:6,1) = rayEnd';
            stnDist = sqrt((refLoc(1,2)-endLoc(1,2))^2 + ...
                (refLoc(1,3)-endLoc(1,3))^2);
            rayAttNow(7,1) = stnDist;
            rayAttNow(8,1) = findSlopeNew(refLoc(1,2:3),endLoc(1,2:3));
            % compute for causal, acusal and symmetric CC
            ccPair = ccStoreFinal(:,j);
            for s = 1:1:length(caseString)
                caseType = caseString(s,1);
                
                plotNum = (s-1)*3+1;
                
                [fExt,distOut,tPU,velPickUpdated(:,s),phPick,spectSNR(:,s),figOut] = ...
                    FTANNewUpd(ccPair,lagVal,minAlpha,maxAlpha,fSamp,fExt,stnDist,...
                    plotSet,caseType,plotNum,Beta);
                %[spectralSNR] = getSpectralSNR(ccPair,fExt,alphaVal(1,filterNo),fSamp,plotNum);
                disp(['Distance is = ', num2str(distOut)]);
            end
            figure(22)
            hold on;
            %plot(fExt,velPickUpdated(:,1),'-o','color','b');
            %plot(fExt,velPickUpdated(:,2),'-o','color','r');
            plot(fExt,velPickUpdated(:,3),'-o','color','g');
            ylim([200,3000])
            hold off;
            % get the observations where v_causal = v_acausal
            goodF = 1;
            for fNo = 1:1:length(fExt)
                velDiff = velPickUpdated(fNo,1)-velPickUpdated(fNo,2);
                velDiffPrct = abs(velDiff/velPickUpdated(fNo,3)*100);
                %% check CC quality
                if(velDiffPrct<20) % acusal and causal velocity should be close
                    if(spectSNR(fNo,3)>3) % check the symmetric SNR
                        if((velPickUpdated(fNo,3)>800)&&(velPickUpdated(fNo,3)<3000))
                            % this obs is good, store it
                            rayStore(ccCount).vel(goodF,1) = velPickUpdated(fNo,3); % stoere the symm velocity
                            rayStore(ccCount).freq(goodF,1) = fExt(fNo,1);
                            rayStore(ccCount).snr(goodF,1) = spectSNR(fNo,3);
                            goodF = goodF +1;
                        end
                    end
                end
            end
            if(goodF>1)
                % at least one freq bin was good
                % hence store the rayAttribute
                rayAttStore(:,ccCount) = rayAttNow;
                ccCount = ccCount+1;
            end
        end
    end
end

% plot the results
% plot the velocities picked
velStore = NaN*ones(length(fExt),ccCount);
snrStore = NaN*ones(length(fExt),ccCount);
for i = 1:1:length(rayStore)
    if(rayAttStore(7,i)>2000)
        for j = 1:1:length(rayStore(i).freq)
            fInd = find(fExt==rayStore(i).freq(j,1));
            velStore(fInd,i) = rayStore(i).vel(j,1);
            snrStore(fInd,i) = rayStore(i).snr(j,1);
        end
    end
end
    
for i = 1:1:length(velStore(1,:))
    figure(1);hold on;
    plot(fExt,velStore(:,i),'b');
    
    %figure(2);hold on;
    %plot(fExt(i,1),snrStore(i,:),'b*');
end