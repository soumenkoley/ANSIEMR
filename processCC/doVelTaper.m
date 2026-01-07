function [newCC,tVec] = doVelTaper(inputCC,fSamp,tVec,vMin,vMax,vInt,dist)
% this function, extends the time axis, dows dense resampling, applies
% velocity taper, and returns back the tapered CC with original sampling
% frequency and time length
dt = 1/fSamp;
tTemp = -100:dt:100;
t1Ind = find(tTemp>=min(tVec),1,'first');
t2Ind = find(tTemp>=max(tVec),1,'first');
newCC = zeros(length(tTemp),1);
newCC(t1Ind:t2Ind) = inputCC;
% now interpolate it;
% resample at 250 Hz
fSampNew = 250;
dtNew = 1/fSampNew;
tNew = -100:dtNew:100;
newCC = interp1(tTemp,newCC,tNew);
% do the velTaper
[w] = vel_taper2(tNew,1/fSampNew,dist/1000,vMin,vMax,vInt);
newCC = newCC'.*(w');
% interp newCC along the old axis
DSFact = fSampNew/fSamp;
newCC = downsample(newCC,DSFact);
newCC = interp1(tTemp,newCC,tVec');
%disp('Done!');
end

