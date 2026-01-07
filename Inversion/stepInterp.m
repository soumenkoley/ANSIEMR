function [yIntp] = stepInterp(x,y,dx,yIntp)
%This function does a step interpolation for a velocity model like
% input, x contains the heights, y the velocities
% xIntp is the final x interval on whcih the interpolation is to be done

yAll = [];
xAll = [];
for i = 2:1:length(x)
    hTemp = x(i-1,1):dx:x(i);
    if(i>=2)
        hTemp = (x(i-1,1)+0.1):dx:x(i);
    end
    nHTemp = length(hTemp);
    yAll = [yAll;y(i,1)*ones(nHTemp,1)];
    xAll = [xAll;hTemp'];
end
yIntp = interp1(xAll,yAll,yIntp);

end

