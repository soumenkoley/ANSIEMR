function [out] = getLineAngle(A,B)
%This function was written to get the anticlockwise angle of a line with respect to
%+ve x-axis
if(B(1,2)>A(1,2))
    if(B(1,1)>A(1,1))
        out = atand((B(1,2)-A(1,2))/(B(1,1)-A(1,1)));
    else
        out = 180 + atand((B(1,2)-A(1,2))/(B(1,1)-A(1,1)));
    end
else
    if(B(1,1)<A(1,1))
        out = 180+atand((B(1,2)-A(1,2))/(B(1,1)-A(1,1)));
    else
        out = 360+atand((B(1,2)-A(1,2))/(B(1,1)-A(1,1)));
    end
end

