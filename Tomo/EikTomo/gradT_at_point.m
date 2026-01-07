function [Tx, Ty] = gradT_at_point(T, x, y, x0, y0, h)
% Compute gradient of T at continuous (x,y) using central differences
% via bilinear interpolation of T itself.

eps = 0.5 * h;   % step for derivative, can be smaller if you like

Tpx = interp_bilinear(T, x+eps, y,   x0, y0, h);
Tmx = interp_bilinear(T, x-eps, y,   x0, y0, h);
Tpy = interp_bilinear(T, x,   y+eps, x0, y0, h);
Tmy = interp_bilinear(T, x,   y-eps, x0, y0, h);

Tx = (Tpx - Tmx) / (2*eps);
Ty = (Tpy - Tmy) / (2*eps);
end