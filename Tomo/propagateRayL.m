function [xNew,yNew,newBox,oldBox,lineLength] = propagateRayL(x,y,currBox,slopeLine,nBoxX,dx,dy,xEnd,yEnd,endBox)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nY = floor(currBox/nBoxX);
nX = rem(currBox,nBoxX);

if(nX == 0)
    currBoxCornerX = (nBoxX-1)*dx;
    currBoxCornerY = (nY)*dy;
else
    currBoxCornerX = (nX-1)*dx;
    currBoxCornerY = (nY+1)*dy;
end

cornerSlope = atand((y-currBoxCornerY)/(x-currBoxCornerX));
if(cornerSlope<0)
    cornerSlope = cornerSlope+180;
end

if(currBox==endBox)
    xNew = xEnd;
    yNew =yEnd;
    newBox = currBox;
    oldBox = currBox;
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
    return;
end

if(slopeLine>cornerSlope)
    newBox = currBox-1;
    oldBox = currBox;
    xNew = currBoxCornerX;
    yNew = y+(-currBoxCornerX+x)*tand(180-slopeLine);
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
elseif(slopeLine==cornerSlope)
    newBox = currBox+nBoxX-1;
    oldBox = currBox;
    yNew = currBoxCornerY;
    xNew = currBoxCornerX;
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
elseif(slopeLine<cornerSlope)
    newBox = currBox+nBoxX;
    oldBox = currBox;
    yNew = currBoxCornerY;
    xNew = x-(currBoxCornerY-y)/tand(180-slopeLine);
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
end
        
end

