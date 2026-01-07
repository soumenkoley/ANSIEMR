% script written to compare star dispersion curves
clear; close all;
starNum = [40,42,44,43,45];
fPathA = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Star';
fPathB = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\NewStars\Star';
colVal= ['b','r','k','m','g','y'];
figure(1);hold on;
for i = 1:1:length(starNum)
    fP = [fPathA,num2str(starNum(i)),'\Fund.txt'];
    A = load(fP);
    plot(A(:,1),A(:,2),'color',colVal(i),'DisplayName',['Star',num2str(starNum(i))]);
end

starNumB = 50;
for i = 1:1:length(starNumB)
    fP = [fPathB,num2str(starNumB(i)),'\Fund.txt'];
    A = load(fP);
    plot(A(:,1),A(:,2),'color',colVal(i),'LineStyle','--',...
        'DisplayName',['Star',num2str(starNumB(i)),'New']);
end

hold off;
legend