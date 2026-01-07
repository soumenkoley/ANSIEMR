function [beamLocs] = plotDoughnut(beamMat,pArray,phiArray,fVec,freqNo)
% get doughnut plot
% this script creates the doughnut plot from the beammatrix
lPArray = length(pArray);
lPhiArray = length(phiArray);
for i = 1:1:lPhiArray
    sI = (i-1)*lPArray+1;
    eI = sI + lPArray-1;
    beamLocs(sI:eI,1) = pArray*cosd(phiArray(i));
    beamLocs(sI:eI,2) = pArray*sind(phiArray(i));
end
beamLocs(:,3) = beamMat(:);
lArrayResponseLocs = length(beamLocs(:,1));
% do something to plot a doughnut accurately enough
% use lenPArray, and lenPhiArray
phiInt = 1;
smallPhi = 10;
startRow = 1;
doughnutInt = (smallPhi/phiInt)*lPArray;
noOfSections = 360/(4*smallPhi/5);
%hold off


figure(10)
subplot(2,3,freqNo)
hold on;
for i = 1:1:(noOfSections)
    rowInt = startRow:1:(startRow+doughnutInt-1);
    tri = delaunay(beamLocs(rowInt,1),beamLocs(rowInt,2));
    
    trisurf(tri,beamLocs(rowInt,1),beamLocs(rowInt,2),beamLocs(rowInt,3));
    shading interp;
    startRow = startRow + doughnutInt-doughnutInt/5;
    if(rowInt(1,end) == lArrayResponseLocs)
        break;
    end
    if((startRow + 2*doughnutInt-doughnutInt/5)>lArrayResponseLocs)
        doughnutInt = lArrayResponseLocs - startRow+1;
    end
end
colorbar;
xlabel('pCos\Phi','FontWeight','bold','FontSize',20);
ylabel('pSin\Phi','FontWeight','bold','FontSize',20);
title([num2str(fVec(freqNo)),' Hz'], 'FontWeight','bold', 'FontSize',20);
set(gca,'FontSize',12);
grid on;
end
