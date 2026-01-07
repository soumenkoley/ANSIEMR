% this script was written to see the inversion output per grid
% it goes to VelGridInvertFeb/gridNo/ and then shows the inversion output

clear;
close all;

gridStart = 549; gridEnd = 575;
gridSee = 557;
fPathInverOut = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvertFeb\grid';
% first plot the grid location both on map, and elevation
load('AllGridAlti.mat');
load('newNodeLocs.mat');
load('AllXYGridPoints.mat');

nX = length(xCoord);
nY = length(yCoord);
% the model matrix will be (nY x nX) matrix
% travel across x, then y to create all set of points

gC = 1;
for i = 1:1:nY
    for j = 1:1:nX
        xyAllGrid(gC,1) = xCoord(j);
        xyAllGrid(gC,2) = yCoord(i);
        gC = gC+1;
    end
end

figure(1);
subplot(1,2,1);
plot(newLoc(:,2),newLoc(:,3),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;
plot(xyAllGrid(gridSee,1),xyAllGrid(gridSee,2),'ko','MarkerSize',8,'MarkerFaceColor','k');
hold off;

subplot(1,2,2);
plot(xyAllGrid(gridStart:gridEnd,1),altIntp(gridStart:gridEnd,1),'bo','MarkerFaceColor','b','MarkerSize',6);
hold on;
plot(xyAllGrid(gridSee,1),altIntp(gridSee,1),'ko','MarkerSize',8,'MarkerFaceColor','k');
hold off;

%%
% load the inversion output
load([fPathInverOut,num2str(gridSee),'\','invOutMarch.mat']);
% add zeros on top for the staircase plot
vsNew = [zeros(1,500);vs];
vpNew = [zeros(1,500);vp];
rhoNew = [zeros(1,500);rho];
hNew = [zeros(1,500);h];

% time to get the mean model
% first we interpolate
hIntp = 0:1:1000;
for i = 1:1:length(vp(1,:))
    vsIntp(:,i) = stepInterp(hNew(:,i),vsNew(:,i),1,hIntp');
    vpIntp(:,i) = stepInterp(hNew(:,i),vpNew(:,i),1,hIntp');
    rhoIntp(:,i) = stepInterp(hNew(:,i),rhoNew(:,i),1,hIntp');
end
% now get mean
vsMean = mean(vsIntp,2);
vpMean = mean(vpIntp,2);
rhoMean = mean(rhoIntp,2);

figure(2);
subplot(1,3,1);
stairs(vpNew,hNew,'b');
hold on;
stairs(vp(:,1),h(:,1),'k','LineWidth',2);
plot(vpMean,hIntp,'r','LineWidth',2);
hold off;
ylabel('Depth (m)');
xlabel('Vp (m/s)');
set(gca,'YDir','reverse');

subplot(1,3,2);
stairs(vsNew,hNew,'b');
hold on;
stairs(vs(:,1),h(:,1),'k','LineWidth',2);
plot(vsMean,hIntp,'r','LineWidth',2);
hold off;
ylabel('Depth (m)');
xlabel('Vs (m/s)');
set(gca,'YDir','reverse');

subplot(1,3,3);
stairs(rhoNew,hNew,'b');
hold on;
stairs(rho(:,1),h(:,1),'k','LineWidth',2);
plot(rhoMean,hIntp,'r','LineWidth',2);
hold off;
ylabel('Depth (m)');
xlabel('Density (kg/m^3)');
set(gca,'YDir','reverse');

%%
% now plot the phase velocities and group velocities for the inversion
fVec = 1.6:0.05:5;
figure(3);
subplot(1,2,1);
plot(fVec,1./vPhTh,'b');
hold on;
plot(fVec,1./vPhTh(:,1),'r','LineWidth',2);
hold off;
xlabel('Frequency (Hz)');
ylabel('Phase velocity (m/s)')

subplot(1,2,2);
plot(fVec,1./vGrpTh,'b');
hold on;
plot(fVec,1./vGrpTh(:,1),'r','LineWidth',2);
hold off;
xlabel('Frequency (Hz)');
ylabel('Group velocity (m/s)')

%%
% plot the relative error
figure(4);
hold on;
plot(relErr,'bo','MarkerSize',6,'MarkerFaceColor','b');
%plot(1:length(relErr),min(relErr),'r--','LineWidth',2);
hold off;