function [vgOut] = getObsGrp(cObs)
%this function converts the phase velocity to group velocity
% the first column in the input is frequency and the second column is phase
% velocity in m/s
%now do u wanna interpolate it and then differentiate
fInt = (1.2:0.1:max(cObs(:,1)))';
phVelInt = interp1(cObs(:,1),cObs(:,2),fInt,'linear','extrap');

% vg = v_phase(1-w/v_phase*dv_phase/dw)^-1; % read Boschi et al 2014

phVelDerivative = gradient(phVelInt,(fInt(2,1)-fInt(1,1)));
vgBC = phVelInt(1:(end),1)./(1-(fInt(1:(end),1)./phVelInt(1:(end),1)).*phVelDerivative);

vgOut = [fInt,vgBC];
%figure(1)
%hold on;
%plot(fInt,phVelInt,'b');
%plot(fInt(1:(end),1),vgBC,'r');

end

