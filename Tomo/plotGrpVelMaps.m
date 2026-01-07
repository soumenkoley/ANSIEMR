% this script was written to load a group velocity map, and see
% which grid points have velocities that might correspond to the higher
% mode

clear; close all;

% load the velocity limits
load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\TheoreticalGrpMaps\thGrpVelOvertLimits.mat');

grpMapPath = 'B:\LimburgBigSurvey1CC-Pair\VelMapsJan\';
fExt = 1.6:0.05:5.0;
fSee = 2;
fSeeInd = find(fExt>=fSee,1,'first');
fName = ['grpVel',num2str(fSee),'.mat'];

load([grpMapPath,fName]);
for i = 1:1:length(xCoord)
    for j = 1:1:length(yCoord)
        vNow = 1/mFinalMatNew(j,i);
        %if(vNow<thGrpVel(fSeeInd,1))
        if(vNow>800)
            mFinalMatNew(j,i) = NaN;
        end
    end
end
figure(1)
h = imagesc(xCoord,yCoord,1./mFinalMatNew);
set(h, 'AlphaData', ~isnan(mFinalMatNew))
shading interp;
colormap('hsv');
colorbar;
xlabel('Distance along Longitude (m)');
ylabel('Distance along Latitude (m)');
titleStr = (['f = ',num2str(fSee),' Hz']);
title(titleStr);
set(gca,'YDir','normal');
hold on;
%plot3(newLoc(:,2),newLoc(:,3),...
%    1000*ones(length(newLoc(:,1)),1),'ko','MarkerSize',3,'MarkerFaceColor','k')