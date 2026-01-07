function [xlong,ylat] = calculatedist(locA,locB,R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
deltaLat = (locA(1,1)-locB(1,1))*pi/180; % in radians converted
deltaLong = (locA(1,2)-locB(1,2))*pi/180;
LatAvg = ((locA(1,1)+locB(1,1))/2)*pi/180;

xlong = R*deltaLong*cos(LatAvg);
ylat = R*deltaLat;

end

