function [AB] = convPh2Grp(f,c)
% this script converts phase to group velocity
%fInt = (1.6:0.1:3.5)';
fInt = f;
A(:,1) = f';
A(:,2) = c;
phVelInt = interp1(A(:,1),A(:,2),fInt,'linear','extrap');
%phVelInt = smooth(phVelInt);
% vg = v_phase(1-w/v_phase*dv_phase/dw)^-1; % read Boschi et al 2014

phVelDerivative = gradient(phVelInt,(fInt(2,1)-fInt(1,1)));
vgBC = phVelInt(1:(end),1)./(1-(fInt(1:(end),1)./phVelInt(1:(end),1)).*phVelDerivative);

AB = [fInt,vgBC];
% figure(100)
% hold on;
% plot(fInt,phVelInt);
% 
% figure(101)
% hold on;
% plot(fInt(1:(end),1),smooth(vgBC));

%disp('Check');
end

