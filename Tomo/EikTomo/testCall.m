clear; close all;

Nx = 40; Ny = 40;
h  = 200;              % 200 m
v0 = 2000;             % m/s
% Low velocity inside (slow anomaly)
v_anom = 500;               % m/s (slower than background)
R   = 1000;                  % radius of anomaly [m]
step_frac = 0.5;
max_steps = 2000;
tol_dist  = 0.5*h;

%%
s  = (1/v0) * ones(Ny, Nx);
x0 = 0; y0 = 0;
x = x0 + (0:Nx-1)*h;
y = y0 + (0:Ny-1)*h;
[X, Y] = meshgrid(x, y);

% Put source somewhere inside the domain, *not* on a node
% xs = 3000;   % m
% ys = 3000;   % m
% xr = 4200;  yr = 1200;    % receiver
xs = x(5);       ys = y(10);             % near top-left-ish
xr = x(end-5);   yr = y(end-10);         % near bottom-right-ish
%%
% add an anomaly to velocity model
% Center and radius of anomaly
xc  = 0.7*(x(1) + x(end));   % center in x (middle of model)
yc  = 0.7*(y(1) + y(end));   % center in y

s_anom = 1/v_anom;

dist_center = sqrt( (X - xc).^2 + (Y - yc).^2 );
mask = dist_center <= R;

s(mask) = s_anom;            % overwrite slowness in anomaly

%% visualize this anomaly
V = 1 ./ s;   % velocity grid

figure(1);
imagesc(x, y, V);
set(gca, 'YDir', 'normal');
axis equal tight;
colorbar;
title('Velocity model (m/s)');
xlabel('x (m)'); ylabel('y (m)');

%%
T = fmm2d_cont_source(s, xs, ys, x0, y0, h);

figure(2)
imagesc(T); axis equal tight; colorbar;
title('Travel time field (constant v, off-node source)');

ray = backtrack_ray_hybrid(T, x0, y0, h, xs, ys, xr, yr, ...
                                step_frac, max_steps, tol_dist);

Lrow = ray_to_lengths(ray, x0, y0, h, Ny, Nx);
%% plot the ray on the contours
figure(3);
contourf(X, Y, T, 20); colorbar; hold on;
plot(xs, ys, 'r*', 'MarkerSize', 10, 'LineWidth', 2);  % source
plot(xr, yr, 'ko', 'MarkerSize', 8, 'LineWidth', 2);   % receiver
plot(ray(:,1), ray(:,2), 'w-', 'LineWidth', 2);        % ray path
%axis equal tight
title('Eikonal travel-time field with backtracked ray');
%% convert the Lrow back to matrix
bb = reshape(Lrow,Ny,Nx);
figure(4)
imagesc(x,y,bb);
