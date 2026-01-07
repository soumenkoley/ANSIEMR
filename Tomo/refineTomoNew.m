% this script performs the tomograpy at any period or frequency

clear;
close all;
% load the rayStore file
%load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\GroupVelPicks\GVPickTrynew2.mat');
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
%load('D:\GVPicksPassive2\AllGVPickUpdt.mat');
fGet = 1.4; % units in Hz
fPathA = 'B:\LimburgBigSurvey1CC-Pair\FinalGrpVelPics\';
load([fPathA,'grpTime',num2str(fGet),'.mat']);
%load('D:\GreensFunction\UpdatedPeriodGroupVelPicks\2_7Hzpick.mat');
% where to store, the ouput maps
fPathStore= 'B:\LimburgBigSurvey1CC-Pair\GrpVelMapsFinal\';
dx = 100; dy = 100;
minRayLength = 1; % units in meters
sigmaVal = 200;
lambdaVal = 0.5;
alpha = 5000;
alphaNew = 4000;
beta = 300;
fExt = fGet;
fGetInd = find(fExt>=fGet,1,'first');

[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);

% specify the name of the file
fStoreName = ['grpVel',num2str(fGet),'.mat'];
figure(1) 
subplot(1,2,1)
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;

rayAttributeStoreFull = rayAttributeStoreFull';
[s1,s2] = size(rayAttributeStoreFull);
nBoxes = (floor(xLimit/dx))*(floor(yLimit/dy));
rayCount = zeros(nBoxes,1);
for rayNo = 1:1:s1
    if(~isnan(phaseVelStoreFull(1,rayNo)))
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
meanPhaseVel = zeros(1,1);
for i = 1:1:1
    notNaNCount = 0;
    for j = 1:1:(l2-l1+1)
        if(phaseVelStoreFull(i,j)~=0)
            if((~isnan(phaseVelStoreFull(i,j))))
                meanPhaseVel(i,1) = 1/phaseVelStoreFull(i,j)+meanPhaseVel(i,1);
                notNaNCount = notNaNCount+1;
            end
        end
    end
    meanPhaseVel(i,1) = meanPhaseVel(i,1)/(notNaNCount);
end
% no need to convert it to slowness
%meanGroupVel = 1./meanGroupVel; % converted it to slowness;
figure(2)
subplot(1,3,1)
plot(fGet,meanPhaseVel','ro','LineWidth',3);
% now compute te kernel matrix G, in d = Gm
[G,rayCountObs,tObs,rayAttUse] = createRayKernel(phaseVelStoreFull,rayAttributeStoreFull,...
                                 timePickStoreFull,dx,dy,fGet,minRayLength);
nPar = nBoxes;
m0 = meanPhaseVel(fGetInd,1)*ones(nPar,1);
%m0 = (1/250)*ones(nPar,1);
tPred = G*m0;
deltaT = tObs-tPred;
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

% compute the velocity distribution
velDistribution = rayAttUse(:,3)./tObs;
figure(6)
subplot(2,1,1)
plot(rayAttUse(:,8),velDistribution,'bo');
subplot(2,1,2)
plot(tObs,velDistribution,'bo');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now here you remove the bad measurements and redo the tomography
stdDeltaTFinal = std(deltaTFinal);
goodCount = 1;
for i = 1:1:length(deltaTFinal)
    if(abs(deltaTFinal(i,1))<6*stdDeltaTFinal)
        deltaTFinalUse(goodCount,1) = deltaT(i,1);
        rayAttUseNew(goodCount,:) = rayAttUse(i,:);
        tObsNew(goodCount,1) = tObs(i,1);
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

% get resolution matrix
resMat = ((GNew'*((CNew^(-1))*GNew)+QNew)^(-1))*(GNew'*(CNew^(-1))*GNew);

for i = 1:1:nPar
    resRow = resMat(i,:);
    subResMat = vec2mat(resRow,nBoxRow);
    [resRadius] = getResolution(subResMat,dx,dy);
    resRadiusVec(1,i) = resRadius;
end
resRadiusMat = vec2mat(resRadiusVec,nBoxRow);
figure(229)
surf(xCoord,yCoord,resRadiusMat);
shading interp;
hold on;
plot3(newLoc(:,2),newLoc(:,3),...
    1000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',6,'MarkerFaceColor','k');

% now plot the output model by converting it back into velocity
mFinalMatNew = vec2mat(mFinalNew,nBoxRow);
timeCorrNew = GNew*mPerturbNew;
deltaTFinalNew  = deltaTFinalUse-GNew*mPerturbNew;
stdFinal = std(deltaTFinalNew);

figure(10)
imagesc(xCoord,yCoord,1./mFinalMatNew);
shading interp;
colormap('jet');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['f = ',num2str(fGet),' Hz, travel time residual = ',num2str(stdFinal),' s']);
title(titleStr);
hold on;
plot3(newLoc(:,2),newLoc(:,3),...
    1000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',6,'MarkerFaceColor','k');
figure(9)
subplot(1,2,1)
imagesc(xCoord,yCoord,1./mFinalMatNew);
shading interp;
colormap('jet');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['f = ',num2str(fGet),' Hz, travel time residual = ',num2str(stdFinal),' s']);
title(titleStr);
set(gca,'YDir','normal');
hold on;
plot3(newLoc(:,2),newLoc(:,3),...
    6000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',6,'MarkerFaceColor','k')
view(2);
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
minVel = min(velVec)-50;
maxVel = max(velVec)+50;

colorVec = jet(300);

figure(9)
subplot(1,2,2)
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'ko','MarkerFaceColor','r','MarkerSize',6);
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
colormap('jet');
caxis([minVel,maxVel]);
hold off;
disp(['RayCount = ', num2str(length(tObsNew))]);

%save([fPathStore,fStoreName],'mFinalNew');