% compare pahse velocity curves
% this script ws written to compare the phase velocity curves of adjacent
% arrays
clear; close all;
fPath = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\';


A = load([fPath,'subArray00\phaseVelNew.mat']);
B = load([fPath,'subArray08\phaseVel.mat']);
C = load([fPath,'subArray05\phaseVelNew.mat']);
D = load([fPath,'subArray13\phaseVel.mat']);
E = load([fPath,'subArray09\phaseVelNew.mat']);
F = load([fPath,'subArray10\phaseVelNew.mat']);
G = load([fPath,'subArray11\phaseVelNew.mat']);

figure(1);
hold on;
plot(A.xi,A.yi,'b');
plot(B.xi,B.yi,'r');
plot(C.xi,C.yi,'g');
plot(D.xi,D.yi,'k--');
plot(E.xi,E.yi,'k');
plot(F.xi,F.yi,'c');
plot(G.xi,G.yi,'m');
