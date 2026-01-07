function Fxy = interp_bilinear(F, x, y, x0, y0, h)
% INTERP_BILINEAR bilinear interpolation of F(i,j) at continuous (x,y)
% Grid: x(j) = x0 + (j-1)*h, y(i) = y0 + (i-1)*h

[Ny, Nx] = size(F);
% fractional indices
jx = (x - x0) / h + 1;
iy = (y - y0) / h + 1;

% clamp to domain
if jx < 1,     jx = 1;     end
if jx > Nx-1,  jx = Nx-1;  end
if iy < 1,     iy = 1;     end
if iy > Ny-1,  iy = Ny-1;  end

j0 = floor(jx); j1 = j0 + 1;
i0 = floor(iy); i1 = i0 + 1;

tx = jx - j0;
ty = iy - i0;

% corner values
F00 = F(i0,j0);
F10 = F(i1,j0);
F01 = F(i0,j1);
F11 = F(i1,j1);

% bilinear interpolation
Fxy = (1-tx)*(1-ty)*F00 + tx*(1-ty)*F01 + (1-tx)*ty*F10 + tx*ty*F11;

end