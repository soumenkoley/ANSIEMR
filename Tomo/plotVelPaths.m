% this script was written to draw virtual velocity paths, which are the
% rays from one particular source to all the receivers and with the
% colorbar showing the apparent velocity of propagation

clear;
close all;
% load the rayStore file
load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\GroupVelPicks\GVPickTryNew5.mat');
load('C:\Dropbox\EinsteinTelescopeSurvey\GreensFunctionPassive2\nodeLocationsCartesian.mat');


sourceNode = 150;
freq = 5; % units in Hz

fExt = 2.6:0.1:8.0;
fInd = find(fExt==freq);

sourceNodeIndA = find(rayAttributeStore(:,1)==sourceNode);
sourceNodeIndB = find(rayAttributeStore(:,2)==sourceNode);

sourceNodeInd = [sourceNodeIndA;sourceNodeIndB];

velStore = groupVelStore(fInd,sourceNodeInd);
maxVel = max(velStore)+5; minVel = min(velStore)-5;

RayStartStore = rayAttributeStore(sourceNodeInd,4:5);
RayEndStore = rayAttributeStore(sourceNodeInd,6:7);

colLevel = 300;
colMat = jet(colLevel);
figure(1)
plot(nodeLocationsCartesian(:,2),nodeLocationsCartesian(:,3),'bo');
hold on;

for i = 1:1:length(sourceNodeInd)
    if(~isnan(velStore(1,i)))
        rayStart = [RayStartStore(i,1);RayEndStore(i,1)];
        rayEnd = [RayStartStore(i,2);RayEndStore(i,2)];
        
        colInd = floor((velStore(1,i)-minVel)/(maxVel-minVel)*colLevel);
        colNow = colMat(colInd,:);
        plot(rayStart,rayEnd,'color',colNow);
    end
end
colorbar;
colormap(jet(colLevel));
caxis([minVel,maxVel]);
hold off;