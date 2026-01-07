function [elements,centers,meanVal,modeVal] = produceHist(dataCol,nBins,lowerLimit,upperLimit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

binSize = (upperLimit-lowerLimit)/nBins;
centers = lowerLimit:binSize:(upperLimit-binSize);
elements = zeros(length(centers),1);
meanVal = 0;
counterVal = 1;
for i = 1:1:length(dataCol)
    if(~isnan(dataCol(i,1)))
        dataCurrInd = floor(((dataCol(i,1)-lowerLimit)/binSize)) + 1;
        elements(dataCurrInd,1) = elements(dataCurrInd,1)+1;
        meanVal = meanVal+dataCol(i,1);
        counterVal = counterVal+1;
    end
end
meanVal = meanVal/(counterVal-1);
modeValInd = find(elements == max(elements));
modeVal =  centers(1,modeValInd(1,1));
figure(1)
stem(centers,elements);
end

