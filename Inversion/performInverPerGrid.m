% this script creates the necessary input files to perform thesurface wave
% inversion
% we use the dinver function from the geopsy package

clear; close all;
gridNo = 643;
NBest = 500;
f1 = 1.6;
ifPlot = 1;
%% several path variables
%path to geopsy excutables
gPath = 'B:\geopsypack-win64-3.4.2\bin\';
% path to parameter file
%fPathParam = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\doInver\myNewParam.param';
fPathParam = 'C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\doInver\myNewParamApril.param';
% working path
wPath = 'B:\LimburgBigSurvey1CC-Pair\VelGridInvertApril\';
totWPath = [wPath,'grid',num2str(gridNo),'\'];

saveFName = 'velMod.model';
saveFPath = [totWPath,saveFName];

savePhDispName = 'RayPhDisp.txt';
savePhDispPath = [totWPath,savePhDispName];

saveGrpDispName = 'RayGrpDisp.txt';
saveGrpDispPath = [totWPath,saveGrpDispName];

% create the phase velocity text file
phFundVel = load('PhFundAvg.txt');
%phFundVel = load('phVelNowSmooth.txt');
fVecFundPh = phFundVel(:,1);

%phOvertVel = load('PhMode1FK.txt');
phOvertVel = load('PhMode1FKHard.txt');

%phOvertVel = load('PhMode1FK0511.txt');
%phOvertVel = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\SubArrayAnalysis\subArray0004vCCBF.txt');

%phOvertVel = load('PhMode1FK11.txt');
%phOvertVel = load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\BeamformOutput\OvertFKSmooth.txt');

fVecOvertPh = phOvertVel(:,1);
%%
% change the directory
cd(totWPath)
currDir = pwd;

% we save slowness instead of velocity
% phVelArray = [fVecFundPh,1./phFundVel(:,2)];
% fId = fopen('phVelFundRayleigh.txt','w');
% fprintf(fId,'%6.2f %6.8f\n',phVelArray');
% fclose(fId);

phVelSmooth = load('phVelNowSmooth.txt');
f1Ind = find(phVelSmooth(:,1)>=f1,1,'first');
fVec = phVelSmooth(f1Ind:end,1);
% we save slowness instead of velocity
phVelArray = [fVec,1./phVelSmooth(f1Ind:end,2)];
fId = fopen('phVelFundRayleigh.txt','w');
fprintf(fId,'%6.2f %6.8f\n',phVelArray');
fclose(fId);

phVelOvertArray = [fVecOvertPh,1./phOvertVel(:,2)];
fId = fopen('phVelOvertRayleigh.txt','w');
fprintf(fId,'%6.2f %6.8f\n',phVelOvertArray');
fclose(fId);

% create the group velocity text file
grpVelSmooth = load('velNowSmooth.txt');
fVec = grpVelSmooth(f1Ind:end,1);
% we save slowness instead of velocity
grpVelArray = [fVec,1./grpVelSmooth(f1Ind:end,2)];
fId = fopen('grpVelRayleigh.txt','w');
fprintf(fId,'%6.2f %6.8f\n',grpVelArray');
fclose(fId);

% copy the parameter file, this stays same
if(exist('myNewParamApril.param'))
    delete myNewParamApril.param;
end
copyfile(fPathParam,currDir);

% create the target file
% delete and recreate if existing
if(exist('my.target'))
    delete my.target;
end
% get the command as a string
% command to add phase velocity
cm1 = ['Type phVelFundRayleigh.txt | ',gPath,'gptarget A -dispersion-rayleigh 0 my.target'];

cm2 = ['Type phVelOvertRayleigh.txt | ',gPath,'gptarget A -dispersion-rayleigh 1 my.target'];

% command to add group velocity
cm3 = ['Type grpVelRayleigh.txt | ',gPath,'gptarget A -dispersion-rayleigh 0 -group my.target'];

stat1 = unix(cm1);
%stat2 = unix(cm2);
stat3 = unix(cm3);

if(exist('my.report'))
    delete my.report;
end
% perform the inversion
cm4 = [gPath,'dinver -i DispersionCurve -optimization -target my.target -param myNewParamApril.param -ns0 50 -ns 30000 -nr 50 -o my.report'];
stat4 = unix(cm4);

modOutFile = ['best',num2str(NBest),'.model'];
cm5 = [gPath,'gpdcreport my.report',' -best ',num2str(NBest),' > ',modOutFile];
stat5 = system(cm5);

delete my.report;
load('grpVel.mat');

cd('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\doInver');
[h,vp,vs,rho,relErr] = readModel([totWPath,modOutFile]);

% write the layer parametrs to a model file to compute the dispersion
% using gpdc

writeModel(h,vp,vs,rho,saveFPath);
% get phase dispersion
cm6 = [gPath,'gpdc -R 1 -step 0.05 -s frequency -min 1.0 -max 5.0 ',saveFPath,' > ',savePhDispPath];
stat6 = system(cm6);

% get group dispersion
cm7 = [gPath,'gpdc -R 1 -step 0.05 -s frequency -group -min 1.0 -max 5.05 ',saveFPath,' > ',saveGrpDispPath];
stat7 = system(cm7);

[fPhTh,vPhTh] = readDispModel([totWPath,'RayPhDisp.txt']);
[fGrpTh,vGrpTh] = readDispModel([totWPath,'RayGrpDisp.txt']);

% save the output velocity models and the misfits
%save([totWPath,'\invOutApril.mat'],'h','vs','vp','rho','relErr','vPhTh','vGrpTh','fPhTh','fGrpTh');

if(ifPlot)
    figure(1);
    hold on;
    hNew = [0;h(1:(end),1)];
    for i = 1:1:length(h(1,:))
        stairs([0;vs(:,i)],[0;h(:,i)],'b');
    end
    stairs([0;vs(:,1)],hNew,'k','LineWidth',2);
    set(gca,'YDir','reverse');
    vsNew = [zeros(1,500);vs];
    vpNew = [zeros(1,500);vp];
    rhoNew = [zeros(1,500);rho];
    hNew = [zeros(1,500);h];

    % time to get the mean model
    % first we interpolate
    hIntp = 0:1:1000;
    for i = 1:1:length(vp(1,:))
        vsIntp(:,i) = stepInterp(hNew(:,i),vsNew(:,i),1,hIntp');
        vpIntp(:,i) = stepInterp(hNew(:,i),vpNew(:,i),1,hIntp');
        rhoIntp(:,i) = stepInterp(hNew(:,i),rhoNew(:,i),1,hIntp');
    end
    % now get mean
    vsMean = mean(vsIntp,2);
    vpMean = mean(vpIntp,2);
    rhoMean = mean(rhoIntp,2);
    
    %meanVs = mean(vs,2);
    stairs(vsMean,hIntp,'r','LineWidth',2);
    xlabel('Shear wave velocity (m/s)');
    ylabel('Depth (m)');
    grid on;
    box on;
    hold off;
    
    figure(2);
    hold on;
    plot(1.0:0.05:5,1./vPhTh,'b');
    plot(phVelArray(:,1),1./phVelArray(:,2),'k','LineWidth',2);
    plot(phVelOvertArray(:,1),1./phVelOvertArray(:,2),'m','LineWidth',2);
    hold off;
    
    figure(3);
    hold on;
    plot(1.0:0.05:5,1./vGrpTh,'b');
    plot(grpVelArray(:,1),1./grpVelArray(:,2),'k','LineWidth',2);
    hold off;
    
%     plot(fGrpTh,1./vGrpTh,'b');
%     plot(grpVelArray(:,1),1./grpVelArray(:,2),'k','LineWidth',2);
%     
end

%% plot the grid location

load('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\nodeLocationsCartesian.mat');

minX = min(nodeLocationsCartesian(:,2))-250;
minY = min(nodeLocationsCartesian(:,3))-250;

newLoc(:,1) = nodeLocationsCartesian(:,1);
newLoc(:,2) = nodeLocationsCartesian(:,2)-minX;
newLoc(:,3) = nodeLocationsCartesian(:,3)-minY;

maxX = max(newLoc(:,2));
maxY = max(newLoc(:,3));

fig1 = openfig('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\assembleCrossCorr\Tomo\GrpTimePickFigs\FundVelApril1_6Hz.fig');
figure(fig1);
hold on;
% for i = 1:1:length(newLoc)
%     plot(newLoc(i,2),newLoc(i,3),'bo');
%     %text(newLoc(i,2),newLoc(i,3),['--',num2str(newLoc(i,1))]);
% end
plot3(xyNow(1,1),xyNow(1,2),10000,'k*','MarkerFaceColor','k','MarkerSize',10);
hold off;


