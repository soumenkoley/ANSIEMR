function [resRadius] = getResolution(resMatrix,dx,dy)
    %this function gets the resolution in m given a resolution matrix by
    %fittinga cone to the peak in the matrxi and finding out the radius of
    %it
    
    [a1,a2] = size(resMatrix);
    xVec = (0:dx:((a2-1))*dx);
    yVec = (0:dx:((a1-1))*dy);
    
    % net get the peak
    [maxPerCol,maxRowInd] = max(resMatrix);
    [maxVal,maxColInd] = max(maxPerCol);
    
    maxInd = [maxRowInd(1,maxColInd),maxColInd,];
    maxCoord = [xVec(maxInd(1,2)),yVec(maxInd(1,1))];
    
%     figure(1)
%     surf(xVec,yVec,resMatrix);
%     shading interp;
%     hold on;
%     plot3(maxCoord(1,1),maxCoord(1,2),2*maxVal,'ko','MarkerSize',8,'MarkerFaceColor','k');
    
    % now try to get the resolution radius
    if(maxVal>=0.01) % this is the cutoff for getting the resolution
        % now get the resolution
        % first search along positive x
        resVal = resMatrix(maxInd(1,1),(maxColInd:a2));
        r1Ind = find(resVal<=(maxVal/10),1,'first');
        r1 = r1Ind*dx;
        
        if(isempty(r1))
            r1 = 2*sqrt((dx^2+dy^2)/2);
        end
        
        % the along negative x
        resVal = fliplr(resMatrix(maxInd(1,1),(1:maxColInd)));
        r2Ind = find(resVal<=(maxVal/10),1,'first');
        r2 = r2Ind*dx;
        if(isempty(r2))
            r2 = 2*sqrt((dx^2+dy^2)/2);
        end
        
        % now along increasing y
        resVal = resMatrix((maxInd(1,1):a1),maxColInd);
        r3Ind = find(resVal<=(maxVal/10),1,'first');
        r3 = r3Ind*dx;
        
        if(isempty(r3))
            r3 = 2*sqrt((dx^2+dy^2)/2);
        end
        
        % now along deecreasing y
        resVal = resMatrix((1:maxInd(1,1)),maxColInd);
        r4Ind = find(resVal<=(maxVal/10),1,'first');
        r4 = r4Ind*dx;
        
        if(isempty(r4))
            r4 = 2*sqrt((dx^2+dy^2)/2);
        end
        
        resRadius = sqrt((r1^2 + r2^2 + r3^2 + r4^2)/4);
        minRes = 2*sqrt((dx^2+dy^2)/2);
        
        if(resRadius<minRes)
            resRadius = minRes;
        end
    else
        resRadius = NaN;
    end
%     disp('Trying baby');
end

