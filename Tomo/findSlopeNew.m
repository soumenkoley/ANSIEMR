function [angle] = findSlopeNew(A,B)
% A is the x,y coordinate of point A
% B is the x,y coordinate of point B
m = (A(1,2)-B(1,2))/(A(1,1)-B(1,1));
angle = atand(m);
if(angle<0)
    angle = 180+angle;
end

end