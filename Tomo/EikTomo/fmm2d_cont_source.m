function T = fmm2d_cont_source(s, xs, ys, x0, y0, h)
% FAST_MARCHING_2D_CONTINUOUS_SOURCE
% Solve |?T| = s using Fast Marching with a source at arbitrary (xs, ys).
%
% Inputs:
%   s   = slowness grid [Ny x Nx] (s = 1 / velocity)
%   xs, ys = source coordinates in same units as x0,y0,h
%   x0, y0 = coordinates of grid node (i=1,j=1)
%   h   = grid spacing (dx = dy = h)
%
% Output:
%   T   = travel time field [Ny x Nx]

[Ny, Nx] = size(s);

% Build coordinate arrays for grid nodes
x = x0 + (0:Nx-1) * h;   % 1 x Nx
y = y0 + (0:Ny-1) * h;   % 1 x Ny
[X, Y] = meshgrid(x, y); % Ny x Nx

% Initialize fields
T      = inf(Ny, Nx);
status = zeros(Ny, Nx);   % 0=far, 1=trial, 2=accepted

% ---- 1) Find surrounding nodes for the continuous source ----
% Find indices such that x(j0) <= xs <= x(j0+1), similarly for y
j0 = find(x <= xs, 1, 'last');
if isempty(j0), j0 = 1; end
if j0 >= Nx, j0 = Nx-1; end    % ensure j0+1 exists

i0 = find(y <= ys, 1, 'last');
if isempty(i0), i0 = 1; end
if i0 >= Ny, i0 = Ny-1; end    % ensure i0+1 exists

% The four surrounding nodes (clipped at edges if needed)
neighbors = [
    i0,   j0;
    i0+1, j0;
    i0,   j0+1;
    i0+1, j0+1
];

% Remove duplicates (can happen near edges)
neighbors = unique(neighbors, 'rows');

% ---- 2) Initialize those nodes as accepted with geometric times ----
for k = 1:size(neighbors,1)
    ii = neighbors(k,1);
    jj = neighbors(k,2);

    dxs = X(ii,jj) - xs;
    dys = Y(ii,jj) - ys;
    dist = hypot(dxs, dys);

    % local velocity and travel time from the exact source to this node
    vloc = 1 / s(ii,jj);
    T(ii,jj) = dist / vloc;

    status(ii,jj) = 2;  % accepted
end

% ---- 3) Initialize neighbors of these accepted nodes as trial ----
for k = 1:size(neighbors,1)
    i = neighbors(k,1);
    j = neighbors(k,2);

    for di = -1:1
        for dj = -1:1
            if abs(di) + abs(dj) ~= 1, continue; end  % 4-neighbors only
            ii = i + di; jj = j + dj;
            if ii>=1 && ii<=Ny && jj>=1 && jj<=Nx
                if status(ii,jj) ~= 2   % not accepted yet
                    T(ii,jj) = local_update_continuous(T, s, status, ii, jj, h);
                    status(ii,jj) = 1;  % trial
                end
            end
        end
    end
end

% ---- 4) Main fast marching loop ----
while true
    trial_idx = find(status == 1);
    if isempty(trial_idx)
        break;   % no more trial points: finished
    end

    % Find trial point with smallest T
    [~, kmin] = min(T(trial_idx));
    p = trial_idx(kmin);
    [i, j] = ind2sub([Ny Nx], p);

    % Accept it
    status(i,j) = 2;

    % Update its 4-neighbors
    for di = -1:1
        for dj = -1:1
            if abs(di) + abs(dj) ~= 1, continue; end
            ii = i + di; jj = j + dj;
            if ii>=1 && ii<=Ny && jj>=1 && jj<=Nx
                if status(ii,jj) ~= 2
                    T(ii,jj) = local_update_continuous(T, s, status, ii, jj, h);
                    status(ii,jj) = 1;  % trial
                end
            end
        end
    end
end

end