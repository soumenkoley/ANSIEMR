function [phaseVelOut,fOut] = calculatePhaseVel(phInp,fPhase,grpVel,fVec,dist)
% this functioon was written to get the different phase velocity branches
% folowing equation 11 in Bensen et al 2007
% phInp is not the unwrapped phase
fOut = min(fVec):0.1:max(fVec);
phUse = interp1(fPhase,phInp,fOut);
grpVelUse = interp1(fVec,grpVel,fOut);
tMax = dist./grpVelUse;

N = 0;
for i = 1:1:length(N)
    for j = 1:1:length(fOut)
        den = phUse(j) - pi/4 + 2*pi*N(i);
        c(j,i) = 1/grpVelUse(j) + ((dist*2*pi*fOut(j))^(-1))*den;
        c(j,i) = 1/c(j,i);
    end
end
phaseVelOut = c;
end

