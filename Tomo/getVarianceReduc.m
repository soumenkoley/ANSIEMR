% this script uses the saved results from refinetomo.m and plots the
% variance reduction as a function of frequency
clear; close all;

fVec = 1.6:0.05:5;

fPath = 'B:\LimburgBigSurvey1CC-Pair\VelMapsFeb\';

for i = 1:1:length(fVec)
    fName = [fPath,'grpVel',num2str(fVec(i)),'.mat'];
    load(fName);
    disp('one done');
    err1 = std(deltaT); err2 = std(deltaTFinalNew);
    varRed(i,1) = (err1-err2)/err1*100;
end
figure(1)
plot(fVec,varRed,'b*');
grid on; box on;
xlabel('Frequency (Hz)');
ylabel('% variance reduction');