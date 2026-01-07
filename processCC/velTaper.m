function [taper] = velTaper(tVec,dt,dist,vMin,vMax,deltaV)
% function to get the velocity taper
fSamp = 1/dt;
tS1 = round((dist/(vMax+deltaV))*fSamp);
tS2 = round((dist/(vMax))*fSamp);
tU1 = round((dist/(vMin))*fSamp);
tU2 = round((dist/(vMin-deltaV))*fSamp);

coeff1 = cos_taper(tS1,tS2,1);
coeff1 = fliplr(coeff1);

coeff2 = cos_taper(tU1,tU2,1);

taper = zeros(1,length(tVec));
taper(1,tS1:tS2) = coeff1;
taper(1,tU1:tU2) = coeff2;
taper(1,tS2:tU1) = 1;

end

