function [velMat,velVec] = getCheckBoard(vMean,vSTD,nX,nY,dBox)
    %this function creates the velocity map as a checkerboard model with
    %Vmean as the mean group velocit at that freq, and then perturbed with
    %+-vSTD
    % nX and nY are the number of boxes along X and Y and dBox is the
    % number of cells over which the velocity is constant
    
    velMat = zeros(nY,nX);
    vMax = vSTD; vMin = -vSTD;
    
    vStart = vMax;
    vStartCh = 1;
    for i = 1:1:nY
        
        vChange = 1;
        % define vstart
        vNow = vStart;
        for j = 1:1:nX
            velMat(i,j) = vNow;
            
            xRem = rem(j,(dBox));
            
            if(xRem==0)
                vChange = -vChange;
            end
            
            if(vChange>0)
                vNow = vStart;
            else
                vNow = -vStart;
            end
        end
        if(rem(i,dBox)==0)
            vStartCh = -vStartCh;
        end
        
        if(vStartCh>0)
            vStart = vMax;
        else
            vStart = vMin;
        end
    end
    
%     % now apply a taper to these
%     % define the taper matrix here
%     nHann = floor(nX/dBox)+1;
%     
%     winVec = [];
%     for i = 1:1:nHann
%         winVal = hann(dBox);
%         winVec = [winVec,winVal'];
%     end
%     
%     for i = 1:1:nY
%         winMat(i,:) = winVec;
%     end
%     winMat = winMat(:,1:nX);
%     velMat = velMat.*winMat;
    
    velMat = velMat+vMean;
    
    velVec = [];
    for i = 1:1:nY
        velVec = [velVec;velMat(i,:)'];
    end
end

