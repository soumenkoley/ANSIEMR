clear all;

load('delaysAll0_1to0_2Hz.mat');
load('Velocity0_1to0_2Hz.mat');

colExtract = 1;

goodCount = 1;
for i = 1:1:length(delaysAll)
    if(~isnan(velAll(i,colExtract)))
        dataSet(goodCount,1) = delaysAll(i,1);
        dataSet(goodCount,2) = delaysAll(i,2);
        dataSet(goodCount,3) = delaysAll(i,3);
        dataSet(goodCount,4) = delaysAll(i,4);
        dataSet(goodCount,5) = delaysAll(i,1)/velAll(i,colExtract);
        dataSet(goodCount,6) = velAll(i,colExtract);
        goodCount = goodCount+1;
    end
end

dX = 300; dY = 300;
[positiveLoc,xExtent,yExtent] = setUpGeometry(dX,dY);
distanceMat = zeros(length(dataSet(:,1)),floor(xExtent/dX)*floor(yExtent/dY));

for i = 1:1:length(dataSet(:,1))
%for i = 1:1:1
    nodeAId = dataSet(i,3); nodeBId = dataSet(i,4);
    nodeAIdInd = find(positiveLoc(:,1)==nodeAId);
    nodeALoc = positiveLoc(nodeAIdInd,2:3);
    nodeBIdInd = find(positiveLoc(:,1)==nodeBId);
    nodeBLoc = positiveLoc(nodeBIdInd,2:3);
    [propDetails] = findTravelBoxes(nodeALoc,nodeBLoc,dX,dY,xExtent,yExtent);
    distanceMat(i,propDetails(:,4)') = dataSet(i,6);
    disp('One done');
end

% check which columns are zero
[outputDistMatrix,goodColumns] = checkZeros(distanceMat);
[velEstimate] = findAvg(outputDistMatrix);

[outputVelMat,xAxis,yAxis] = plotVelocityMaps(velEstimate,goodColumns,10,11,dX,dY);

figure(1000)
surf(xAxis,yAxis,outputVelMat);
shading interp;
% % set up the time matrix
% timeMatrix = dataSet(:,5);
% 
% leastSqOperator = outputDistMatrix*outputDistMatrix';
% invOperator  = inv(leastSqOperator);
% slownessEstimate = (outputDistMatrix'*invOperator)*timeMatrix;
% velEstimate = 1./slownessEstimate;

