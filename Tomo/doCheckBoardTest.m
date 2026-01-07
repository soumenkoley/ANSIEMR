% this script performs the tomograpy at any period or frequency

clear;
close all;
% load the rayStore file
%load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\GroupVelPicks\GVPickTrynew2.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
%load('D:\GVPicksPassive2\AllGVPickUpdt.mat');
%load('phaseVelPicks.mat');
%load('D:\GreensFunction\UpdatedPeriodGroupVelPicks\2_7Hzpick.mat');

dx = 200; dy = 200;
minRayLength = 1; % units in meters
sigmaVal = 200;
lambdaVal = 0.2;
alpha = 3000;
alphaNew = 500;
beta = 300;
fGet = 2; % units in Hz
fExt = fGet;
fGetInd = find(fExt >= fGet,1,'first');
load(['C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TravelTimePicks\LowVelAug\SNR7\grpTimeNew',num2str(fGet),'.mat']);[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);

figure(1) 
subplot(1,2,1)
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;

rayAttributeStoreFull = rayAttributeStoreFull';
[s1,s2] = size(rayAttributeStoreFull);
nBoxes = (floor(xLimit/dx))*(floor(yLimit/dy));
rayCount = zeros(nBoxes,1);
for rayNo = 1:1:s1
    if(~isnan(phaseVelStoreFull(fGetInd,rayNo)))
        nodeA = rayAttributeStoreFull(rayNo,1); nodeB = rayAttributeStoreFull(rayNo,2);
        nodeAInd = find(newLoc(:,1)==nodeA); nodeBInd = find(newLoc(:,1)==nodeB);
        nodeALoc = newLoc(nodeAInd,2:3); nodeBLoc = newLoc(nodeBInd,2:3);
        rayStart = [nodeALoc(1,1),nodeBLoc(1,1)]';
        rayEnd = [nodeALoc(1,2),nodeBLoc(1,2)]';
        plot(rayStart,rayEnd,'k','LineWidth',1);
        [propDetails] = findTravelBoxes(nodeALoc,nodeBLoc,dx,dy,xLimit,yLimit);
        rayCount(propDetails(:,4),1) = rayCount(propDetails(:,4),1)+1;
        %plot(propDetails(:,1),propDetails(:,2),'ro','MarkerSize',4,'MarkerFaceColor','r');
    end
end
hold off;

% now get the coordinates of the boxes
boxCoord = zeros(nBoxes,2);
for boxNo = 1:1:nBoxes
    nBoxRow = floor(xLimit/dx);
    nBoxCol = floor(yLimit/dy);
    if(rem(boxNo,nBoxRow)==0)
        colNo = nBoxRow;
        rowNo = floor(boxNo/nBoxRow);
    else
        colNo = rem(boxNo,nBoxRow);
        rowNo = floor(boxNo/nBoxRow)+1;
    end
    boxCoord(boxNo,1) = (colNo-1)*dx + dx/2;
    boxCoord(boxNo,2) = (rowNo-1)*dy + dy/2;
end

rayCountMat = vec2mat(rayCount',nBoxRow);
xCoord = ((1:1:nBoxRow)-1)*dx + dx/2;
yCoord = ((1:1:nBoxCol)-1)*dy + dy/2;

figure(1)
subplot(1,2,2)
surf(xCoord,yCoord,rayCountMat);
%caxis([0,10]);
colormap('jet');
shading interp;

% now compute the mean group velocity
l1 = 1;
l2 = length(phaseVelStoreFull(1,:));

meanPhaseVel = 0;
notNaNCount = 0;
for j = 1:1:(l2-l1+1)
    if(phaseVelStoreFull(fGetInd,j)~=0)
        if((~isnan(phaseVelStoreFull(fGetInd,j))))
            meanPhaseVel = phaseVelStoreFull(fGetInd,j)+meanPhaseVel;
            notNaNCount = notNaNCount+1;
            velTemp(notNaNCount,1) = phaseVelStoreFull(fGetInd,j);
        end
    end
end
%meanPhaseVel = meanPhaseVel/(notNaNCount-1);
%stdVel = 2*std(velTemp);
meanPhaseVel = 250;
stdVel = 100;

disp(['Std = ',num2str(stdVel),' and mean = ',num2str(meanPhaseVel)]);
% no need to convert it to slowness
%meanGroupVel = 1./meanGroupVel; % converted it to slowness;
figure(2)
subplot(1,3,1)
plot(fExt,meanPhaseVel','ro','LineWidth',3);
% now compute the kernel matrix G, in d = Gm
[G,rayCountObs,tObs,rayAttUse] = createRayKernel(phaseVelStoreFull,rayAttributeStoreFull,...
                                 timePickStoreFull,dx,dy,fGet,minRayLength);
% create the theoretical checkerboard model
[vChkMat,vChkVec] = getCheckBoard(meanPhaseVel,stdVel,nBoxRow,nBoxCol,3);
[vChkMatTap,vChkVecTap] = getCheckBoardTapered(meanPhaseVel,stdVel,nBoxRow,nBoxCol,3);

%now convert velocities to slowness
nPar = nBoxes;
m0 = (1/meanPhaseVel)*ones(nPar,1);
slowTh = 1./vChkVec;
tTh = G*slowTh;
tPred = G*m0;
deltaT = tTh-tPred;
stdInitial = std(deltaT);
disp(['Initial rms error = ', num2str(stdInitial)]);
nObs = 1:1:length(tObs);
figure(2)
subplot(1,3,2)
plot(nObs,deltaT,'b*');

subplot(1,3,3)
hold on;
plot(nObs,tObs,'b');
plot(nObs,tPred,'r');

% now get the smoothing matrix
[F] = getSmoothingMat(boxCoord,sigmaVal);
figure(3)
FMat = vec2mat(F(:,34),nBoxRow);
surf(xCoord,yCoord,FMat);
shading interp;

% now get weighting matrix
[H] = getWeightingMat(rayCount,lambdaVal);

% now contruct the covariance matrix
% the error in the observed travel times should be ideally done over
% monthly observations, but here we do something else, set it as deltaT/5
C = abs(diag(deltaT/0.2));
% goodRay = 1;
% for rayNo = 1:1:length(errorTStoreNew(1,:))
%     if(~isnan(errorTStoreNew(fGetInd,rayNo)))
%         sigmaVec(goodRay,1) = errorTStoreNew(fGetInd,rayNo);
%         goodRay = goodRay+1;
%     end
% end
%C = abs(diag(sigmaVec.^2));
% compute the weighting and smoothing matrix equivalent
Q = ((alpha^2)*(F'*F)) + ((beta^2)*(H'*H)) ;
mPerturb = ((G'*((C^(-1))*G)+Q)^(-1))*(G'*(C^(-1))*deltaT);
mFinal = m0+mPerturb;

% now plot the output model by converting it back into velocity
mFinalMat = vec2mat(mFinal,nBoxRow);
figure(4)
surf(xCoord,yCoord,1./mFinalMat);
shading interp;
colormap('jet');
colorbar;
%zlim([50,600]);
%caxis([50,600]);

% lastly compute the error in traveltimes also after inversion
timeCorr = G*mPerturb;
deltaTFinal  = deltaT-G*mPerturb;
figure(5)
subplot(3,1,1)
plot(deltaT);
hold on;
%subplot(2,1,2)
plot(deltaTFinal,'r');
ylim([-0.6,0.6]);

subplot(3,1,2)
plot(timeCorr);

subplot(3,1,3)
plot(1./mFinal,'bo');
stdFinal = std(deltaTFinal);
disp(['Initial rms error = ', num2str(stdFinal)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now here you remove the bad measurements and redo the tomography
stdDeltaTFinal = std(deltaTFinal);
goodCount = 1;
for i = 1:1:length(deltaTFinal)
    if(abs(deltaTFinal(i,1))<3*stdDeltaTFinal)
        deltaTFinalUse(goodCount,1) = deltaT(i,1);
        rayAttUseNew(goodCount,:) = rayAttUse(i,:);
        tObsNew(goodCount,1) = tTh(i,1);
        GNew(goodCount,:) = G(i,:);
        %sigmaVecNew(goodCount,1) = sigmaVec(i,1);
        goodCount = goodCount+1;        
    end
end

% recompute the weighted norm penalizing matrix;
% for that you need to redo the rayCount part
for i = 1:1:length(GNew(1,:))
    rayCountNew(i,1) = length(find(GNew(:,i)>minRayLength));
end

[HNew] = getWeightingMat(rayCountNew,lambdaVal);

% no need to recompute F matrix or the smoothing matrix
% but you need to compute the covariance matrix again

%CNew = abs(diag(sigmaVecNew.^2));
%CNew = diag(0.01*ones(length(deltaTFinalUse),1));
CNew = abs(diag(deltaTFinalUse/0.5));
QNew = ((alphaNew^2)*(F'*F)) + ((beta^2)*(HNew'*HNew)) ;
mPerturbNew = ((GNew'*((CNew^(-1))*GNew)+QNew)^(-1))*(GNew'*(CNew^(-1))*deltaTFinalUse);
mFinalNew = m0+mPerturbNew;

resMat = ((GNew'*((CNew^(-1))*GNew)+QNew)^(-1))*(GNew'*(CNew^(-1))*GNew);
% now handle the resolution mat

% see the 120th node point
for i = 420:1:420
    resRow = resMat(i,:);
    subResMat = vec2mat(resRow,nBoxRow);
end
figure(229)
surf(xCoord,yCoord,subResMat);
shading interp;

% now plot the output model by converting it back into velocity
mFinalMatNew = vec2mat(mFinalNew,nBoxRow);
timeCorrNew = GNew*mPerturbNew;
deltaTFinalNew  = deltaTFinalUse-GNew*mPerturbNew;
stdFinal = std(deltaTFinalNew);

maxVel = meanPhaseVel+stdVel;
minVel = meanPhaseVel-stdVel;

figure(9)
subplot(1,3,1)
imagesc(xCoord,yCoord,vChkMat);
shading interp;
colorbar;
colormap('parula');
caxis([minVel,maxVel]);
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['V_{mean} = ', num2str(meanPhaseVel),' m/s, \sigma = ',num2str(stdVel),' m/s']);
title(titleStr);
hold on;
plot3(newLoc(:,2),newLoc(:,3),...
    1000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',6,'MarkerFaceColor','k')

figure(9)
subplot(1,3,2)
imagesc(xCoord,yCoord,1./mFinalMatNew);
shading interp;
colorbar;
colormap('parula');
caxis([minVel,maxVel]);
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['f = ',num2str(fGet),' Hz, travel time residual = ',num2str(stdFinal),' s']);
title(titleStr);
hold on;
plot3(newLoc(:,2),newLoc(:,3),...
    1000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',6,'MarkerFaceColor','k')
%zlim([50,600]);
%caxis([50,600]);

% lastly compute the error in traveltimes also after inversion

figure(8)
subplot(3,1,1)
plot(deltaTFinalUse);
hold on;
%subplot(2,1,2)
plot(deltaTFinalNew,'r');
ylim([-0.6,0.6]);

subplot(3,1,2)
plot(timeCorrNew);

subplot(3,1,3)
plot(1./mFinal,'bo');
disp(['Final rms error = ', num2str(stdFinal)]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now plot the rayPaths as well with color code varying as avg velocity
% along the path

velVec = rayAttUseNew(:,3)./tObsNew;
minVel = min(velVec)-5;
maxVel = max(velVec)+5;

colorVec = jet(300);

figure(9)
subplot(1,3,3)
hold on;
for i = 1:1:length(tObsNew)
    colorNow = colorVec(floor(((velVec(i,1)-minVel)/(maxVel-minVel))*300),:);
    rayStart = [rayAttUseNew(i,4),rayAttUseNew(i,6)]';
    rayEnd = [rayAttUseNew(i,5),rayAttUseNew(i,7)]';
    plot(rayStart,rayEnd,'color',colorNow);
end
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
title(['Ray Count = ',num2str(length(tObsNew)),' at f = ',num2str(fGet),' Hz']);
colorbar;
colormap('parula');
caxis([minVel,maxVel]);
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ko','MarkerFaceColor','k','MarkerSize',6);
hold off;

%% prepare plot for reduced node list
 load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\allStn.mat');
 notPlotList = ["WKQQA";"W3QZA";"W2RAA";"XFRFA";"W4QGA";"XPQJA";...
    "W7PSA";"W4L7A";"XAL4A";"XELYA";"W9MRA";"XFMHA";"XMMBA";...
    "XFMWA"];
% add plot or not plot thing in the nodeNames cellarray

pCount = 1;
for i = 1:1:length(allStn)
    ss = find(notPlotList == allStn{i,1});
    if(isempty(ss))
        newLocPlot(pCount,1:2) = newLoc(i,2:3);
        pCount = pCount+1;
    end
end
figure(384); hold on;
h = imagesc(xCoord,yCoord,1./mFinalMatNew);
set(h, 'AlphaData', ~isnan(mFinalMatNew))
shading interp;
colormap('parula');
colorbar;
caxis([minVel,maxVel]);
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['f = ',num2str(fGet),' Hz, travel time residual = ',num2str(stdFinal),' s']);
title(titleStr);
set(gca,'YDir','normal');
hold on;
plot3(newLocPlot(:,1),newLocPlot(:,2),...
    1000*ones(length(newLocPlot(:,1)),1),'ko','MarkerSize',3,'MarkerFaceColor','k')
grid on; box on;