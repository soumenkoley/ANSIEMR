function [xNew,yNew,newBox,oldBox,lineLength] = propagateRayR(x,y,currBox,slopeLine,nBoxX,dx,dy,xEnd,yEnd,endBox)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nY = floor(currBox/nBoxX);
nX = rem(currBox,nBoxX);

if(nX == 0)
    currBoxCornerX = nBoxX*dx;
    currBoxCornerY = nY*dy;
else
    currBoxCornerX = nX*dx;
    currBoxCornerY = (nY+1)*dy;
end

cornerSlope = atand((y-currBoxCornerY)/(x-currBoxCornerX));

if(currBox==endBox)
    xNew = xEnd;
    yNew =yEnd;
    newBox = currBox;
    oldBox = currBox;
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
    return;
end

if(slopeLine>cornerSlope)
    newBox = currBox+nBoxX;
    oldBox = currBox;
    yNew = currBoxCornerY;
    xNew = x+(currBoxCornerY-y)/tand(slopeLine);
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
elseif(slopeLine==cornerSlope)
    newBox = currBox+nBoxX+1;
    oldBox = currBox;
    yNew = currBoxCornerY;
    xNew = currBoxCornerX;
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
elseif(slopeLine<cornerSlope)
    newBox = currBox+1;
    oldBox = currBox;
    xNew = currBoxCornerX;
    yNew = y+(currBoxCornerX-x)*tand(slopeLine);
    lineLength = sqrt((x-xNew)^2+(y-yNew)^2);
end
        
end

