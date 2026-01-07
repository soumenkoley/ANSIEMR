% this script performs the tomograpy at any period or frequency

clear;
close all;
% load the rayStore file
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');
%load('D:\GVPicksPassive2\AllGVPickUpdt.mat');
load('phaseVelPicks.mat');
%load('D:\GreensFunction\UpdatedPeriodGroupVelPicks\2_7Hzpick.mat');
dx = 200; dy = 200;
minRayLength = 1; % units in meters
sigmaVal = 200;
lambdaVal = 0.2;
alpha = 500;
beta = 200;
fGet = 1.6; % units in Hz
fExt = 1.2:0.1:2.5;
fGetInd = find(fExt >= fGet,1,'first');
[newLoc,xLimit,yLimit] = setUpGeometry(dx,dy);

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
caxis([0,10]);
colormap('jet');
shading interp;

% now compute the mean group velocity
l1 = 1;
l2 = length(phaseVelStoreFull(:,1));
meanPhaseVel = zeros(length(fExt),1);
for i = 1:1:length(fExt)
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
plot(fExt,meanPhaseVel','ro','LineWidth',3);
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
C = abs(diag(deltaT/1));

% compute the weighting and smoothing matrix equivalent
Q = ((alpha^2)*(F'*F)) + ((beta^2)*(H'*H)) ;
mPerturb = ((G'*((C^(-1))*G)+Q)^(-1))*(G'*deltaT);
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