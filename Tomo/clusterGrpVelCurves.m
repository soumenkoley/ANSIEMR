% this script was written to cluster the group velocity maps in a try to
% separate the two modes out
clear; close all;

fPath = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvert\';

gridNo = 1:1:1122;

gC = 1;

figure(1);
hold on;
for i = 1:1:length(gridNo)
    fName = ['grid',num2str(i)];
    totPath = [fPath,fName,'\grpVel.mat'];
    if(exist(totPath))
        load(totPath);
        if(~isnan(sum(velNow)))
            vAllGrid(:,gC) = smooth(velNow);
            plot(fExt,vAllGrid(:,gC))
            gC = gC+1;
            
        end
    end
end
hold off;



