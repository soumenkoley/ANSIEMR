% this script was written to find out the list of stations lying within a
% given region
clear; close all;

%xv = [50.7767907; 50.7712595; 50.7688168; 50.773485; 50.7837967;50.7823263;50.7767907];
%yv = [5.9313683; 5.9342865; 5.9468355; 5.9549036; 5.9470072; 5.9354023; 5.9313683];

xv = [50.7688168;50.7567591;50.7650116;50.7757051;50.773485;50.7688168];
yv = [5.9468355;5.9632114;5.977631;5.9753136;5.9549036;5.9468355];

refLong = 5.9060966; refLat = 50.7513927;
refPoint = [refLat,refLong];
% convert the edge of the polygon to x,y from latlong
R = 6371000; % units in meters
for i = 1:1:length(xv)
    pointNow = [xv(i,1),yv(i,1)];
    [xvC(i,1),yvC(i,1)] = calculatedist(pointNow,refPoint,R);
end

% columns 1 and 2 of A have the latlong corresponding the station names in allStn 
A = xlsread('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUseNew.xlsx');
% this is string array with the name of the station
[allStn] = loadSensorNames('C:\Users\soume\Dropbox\EinsteinTelescopeSurvey\TerzietBigSurvey1\MetaDataUse.txt');

% convert latlong values in A to carteisan usin the same ref latlong
for i = 1:1:length(A(:,1))
    pointNow = [A(i,1),A(i,2)];
    [xq(i,1),yq(i,1)] = calculatedist(pointNow,refPoint,R);
end

[inInd,onInd] = inpolygon(xq,yq,xvC,yvC);

figure(1)
plot(xq,yq,'bo','MarkerFaceColor','b','MarkerSize',8);
hold on;
plot(xq(inInd,1),yq(inInd,1),'ro','MarkerFaceColor','r','MarkerSize',8);

stnStore = allStn(inInd,1);