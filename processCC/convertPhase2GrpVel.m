% this program computes the group velocity from phase velocity
% load the file first
clear;
%close all;
%fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\';
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\';
%A = load([fPath,'subArray03\','Fund.txt']);
A = load([fPath,'Overt.txt']);
%now do u wanna interpolate it and then differentiate
fInt = (1.2:0.1:2.5)';
phVelInt = interp1(A(:,1),A(:,2),fInt,'linear','extrap');
phVelInt = smooth(phVelInt);
% vg = v_phase(1-w/v_phase*dv_phase/dw)^-1; % read Boschi et al 2014

phVelDerivative = gradient(phVelInt,(fInt(2,1)-fInt(1,1)));
vgBC = phVelInt(1:(end),1)./(1-(fInt(1:(end),1)./phVelInt(1:(end),1)).*phVelDerivative);

AB = [fInt,vgBC];
figure(100)
hold on;
plot(fInt,phVelInt);

figure(101)
hold on;
plot(fInt(1:(end),1),smooth(vgBC));






