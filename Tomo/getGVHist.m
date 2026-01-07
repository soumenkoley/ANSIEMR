% make group velocity histograms
clear; close all;
load('D:\GVPicksPassive2\AllGVPickUpdt.mat');

[s1,s2] = size(groupVelStoreNew);

fExt = 2.6:0.1:8.0;
dV = 5; % in m/s
maxV = 500;
NDV = floor(maxV/dV)+1;
groupVelHist = zeros(s1,NDV);
vVec = (1:1:NDV)*dV;
addGrpVel = zeros(s1,1);
grpVelCount = zeros(s1,1);
for i = 1:1:s1
    notNaNCount = 1;
    for j = 1:1:s2
        if(~isnan(groupVelStoreNew(i,j)))
            indX = floor(groupVelStoreNew(i,j)/dV)+1;
            groupVelHist(i,indX) = groupVelHist(i,indX)+1;
            addGrpVel(i,1) = addGrpVel(i,1)+groupVelStoreNew(i,j);
            tempVel(1,notNaNCount) = groupVelStoreNew(i,j);
            grpVelCount(i,1) = grpVelCount(i,1)+1;
            notNaNCount = notNaNCount+1;
        end
    end
    stdGV(i,1) = std(tempVel);
end

meanGrpVel = addGrpVel./grpVelCount;
figure(1)
imagesc(fExt,vVec,groupVelHist');
shading interp;
hold on;
plot3(fExt,meanGrpVel,1000*ones(length(fExt),1),'k','LineWidth',3);
plot3(fExt,(meanGrpVel+stdGV),1000*ones(length(fExt),1),'k--','LineWidth',3);
plot3(fExt,(meanGrpVel-stdGV),1000*ones(length(fExt),1),'k--','LineWidth',3);
hold off;