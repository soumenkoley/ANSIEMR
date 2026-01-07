function [filterVal] = getGaussFilterFull(f0,alpha,fVec,beta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
filterVal = exp((-alpha*((fVec-f0)/f0).^2));
bandWidth = sqrt(beta/alpha);
omegaHigh = 2*pi*f0*(1+bandWidth);
omegaLow = 2*pi*f0*(1-bandWidth);
omegaVec = 2*pi*fVec;
filterVal((omegaVec<omegaLow)) = 0;
filterVal((omegaVec>omegaHigh)) = 0;
% figure(199)
% fPlot = [1,2,3,4,5,6,7,8];
% fNowInd = find(fPlot==f0);
% 
% if(~isempty(fNowInd))
% hold on;
% plot(fVec,filterVal);
% disp('Filter prepared');
% end
end
