% thi script was written to display he velocity model as a 2D
% along a line passing through some grid points

clear; %close all;

gridStart = 281; gridEnd = 304;

dx = 200; dz = 1; % units in meters
dzMax = 900;
dzArray = 0:dz:dzMax;

gridVal = gridStart:1:gridEnd;
nGrid = length(gridVal);
dxArray = 0:dx:((nGrid-1)*dx);

gPath = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvertFeb\grid';
for i = 1:1:nGrid
    % load the inversion output
    fPathFull = [gPath,num2str(gridVal(i)),'\invOutMarch.mat'];
    if(exist(fPathFull))
        load([gPath,num2str(gridVal(i)),'\invOutMarch.mat']);
        vsMean = mean(vs,2);
        for j = 1:1:length(dzArray)
            zInd = find(h(:,i)>=dzArray(j),1,'first');
            velMat(j,i) = vs(zInd,1);
            velMatMean(j,i) = vsMean(zInd,1);
        end
        %     p = polyfit(dzArray',velMat(:,i),3);
        %     velMatPolyFit(:,i) = p(1)*dzArray.^3 + p(2)*dzArray.^2+ p(3)*dzArray + p(4);
        %     disp('one grid');
        relErrAll(i,1) = min(relErr);
    else
        velMat(:,i) = velMat(:,i-1);
        velMatMean(:,i) = velMatMean(:,i-1);
        relErrAll(i,1) = relErrAll(i-1,1);
    end
end

for i = 1:1:length(dzArray)
    velMat(i,:) = smooth(velMat(i,:),5);
    velMatMean(i,:) = smooth(velMatMean(i,:),5);
end


figure(1);
subplot(1,2,1)
surf(dxArray,dzArray,velMat);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');

figure(1);
subplot(1,2,2)
surf(dxArray,dzArray,velMatMean);
shading interp;
view(2);
set(gca,'YDir','reverse');
colorbar;colormap('hsv');

figure(2);
plot(dxArray,relErrAll);