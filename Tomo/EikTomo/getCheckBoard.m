function [s_true,x,y] = getCheckBoard(vBG,L,amp,Nx,Ny,h)
%checkerboard parameters
%Lx = Ly = L% checker cell size in x (m)
%amp = 0.10;   % ±10% velocity perturbation
% vBg = background velocity (Rayleigh c or S-wave Vs proxy)
% Nx = number of grid points along x
% Ny = number of grid points along y
% h = cell size in x,y
%%
x0 = 0; y0 = 0;
x  = x0 + (0:Nx-1)*h;
y  = y0 + (0:Ny-1)*h;
[X, Y] = meshgrid(x, y);

s_bg = 1 / vBG;      % background slowness
s_true = s_bg * ones(Ny, Nx);

% simple “chessboard” using floor() of positions
ix = floor(X / L);
iy = floor(Y / L);
pattern = (-1).^(ix + iy);     % +1 / -1 alternating pattern

% define velocity perturbation
v_true = vBG .* (1 + amp * pattern);  % ± amp in velocity
s_true = 1 ./ v_true;
end

