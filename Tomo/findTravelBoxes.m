function [propDetails] = findTravelBoxes(nodeA,nodeB,dx,dy,X,Y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nBoxX = floor(X/dx); nBoxY = floor(Y/dy);
% first decide which is the base node
%x1 = nodeA(1,1); y1 = nodeA(1,2); x2 = nodeB(1,1); y2 = nodeB(1,2);
lineSlope = findSlopeNew(nodeA,nodeB);

dARight = X-nodeA(1,1);dBRight = X-nodeB(1,1);

if(lineSlope<=90)
    if(dARight>=dBRight)
        % nodeA and nodeB are in order
        nodeStart = nodeA;
        nodeEnd = nodeB;
    else
        nodeStart = nodeB;
        nodeEnd = nodeA;
    end
else
    if(dARight<=dBRight)
        nodeStart = nodeA;
        nodeEnd = nodeB;
    else
        nodeStart = nodeB;
        nodeEnd = nodeA;
    end
end

x1 = nodeStart(1,1); y1 = nodeStart(1,2); x2 = nodeEnd(1,1); y2 = nodeEnd(1,2);

sBoxA = nBoxX*floor(y1/dy)+floor(x1/dx)+1;
eBoxB = nBoxX*floor(y2/dy)+floor(x2/dx)+1;

if(lineSlope<=90)
    xC = x1; yC = y1;
    currBox = sBoxA;
    oldB = currBox;
    boxCounter = 1;
    while(oldB~=eBoxB)
        [xN,yN,newB,oldB,lineL] = propagateRayR(xC,yC,currBox,lineSlope,nBoxX,dx,dy,x2,y2,eBoxB);
        propDetails(boxCounter,:) = [xN,yN,newB,oldB,lineL];
        xC = xN; yC = yN; currBox = newB;
        boxCounter = boxCounter+1;
    end
end

if(lineSlope>90)
    xC = x1; yC = y1;
    currBox = sBoxA;
    oldB = currBox;
    boxCounter = 1;
    while(oldB~=eBoxB)
        [xN,yN,newB,oldB,lineL] = propagateRayL(xC,yC,currBox,lineSlope,nBoxX,dx,dy,x2,y2,eBoxB);
        propDetails(boxCounter,:) = [xN,yN,newB,oldB,lineL];
        xC = xN; yC = yN; currBox = newB;
        boxCounter = boxCounter+1;
    end
end

if(boxCounter ==1)
    distTravel = sqrt((nodeA(1,1)-nodeB(1,1))^2 + (nodeA(1,2)-nodeB(1,2))^2);
    propDetails(boxCounter,:) = [xC,yC,oldB,oldB,distTravel];
end
%disp('Testing!');

end

