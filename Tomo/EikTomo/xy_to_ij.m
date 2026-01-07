function [i,j] = xy_to_ij(x,y,x0,y0,h,Ny,Nx)
j = max(min(round((x-x0)/h) + 1, Nx),1);
i = max(min(round((y-y0)/h) + 1, Ny),1);
end