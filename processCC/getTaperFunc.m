function [H] = getTaperFunc(tVec,dist,velHigh,velLow)
    %this function gives you the taper function to be applied while
    %generating the theoretical cross-correlation
    
    % following yao and van der hilst 2009
    nT = length(tVec);
    H = zeros(nT,1);
    
    %tLower = delayStn-TW/2; tUpper = delayStn+TW/2;
    tLower = dist/velHigh; tUpper = dist/velLow;
    tLowerInd = find(tVec<=tLower,1,'last');
    tUpperInd = find(tVec>=tUpper,1,'first');
    
    for i = tLowerInd:1:tUpperInd
        phaseNow = 2*pi*(tVec(i,1)-delayStn)/TW;
        H(i,1) = 0.5*(1+cos(phaseNow));
    end
end