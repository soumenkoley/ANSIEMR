function Tnew = local_update_continuous(T, s, status, i, j, h)
% LOCAL_UPDATE_CONTINUOUS Upwind update for isotropic grid (dx = dy = h)
% Uses ONLY accepted neighbors (status==2)

[Ny, Nx] = size(T);

Tx = inf;
Ty = inf;

% X-direction neighbors (accepted only)
if i > 1   && status(i-1,j) == 2, Tx = min(Tx, T(i-1,j)); end
if i < Ny  && status(i+1,j) == 2, Tx = min(Tx, T(i+1,j)); end

% Y-direction neighbors (accepted only)
if j > 1   && status(i,j-1) == 2, Ty = min(Ty, T(i,j-1)); end
if j < Nx  && status(i,j+1) == 2, Ty = min(Ty, T(i,j+1)); end

% If no accepted neighbors, do not update
if isinf(Tx) && isinf(Ty)
    Tnew = T(i,j);
    return;
end

v = 1 / s(i,j);   % local velocity

% If only one accepted neighbor: simple 1D update
if isinf(Tx)
    Tnew = Ty + h / v;
    return;
elseif isinf(Ty)
    Tnew = Tx + h / v;
    return;
end

% Two accepted neighbors: Sethian's quadratic update
a = min(Tx, Ty);
b = max(Tx, Ty);

if (b - a) >= h / v
    % Characteristic cone condition
    Tnew = a + h / v;
else
    % Both neighbors influence the solution
    Tnew = (a + b + sqrt(2*(h/v)^2 - (b - a)^2)) / 2;
end

end