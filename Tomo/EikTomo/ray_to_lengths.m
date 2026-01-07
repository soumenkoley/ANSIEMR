function Lrow = ray_to_lengths(ray, x0, y0, h, Ny, Nx)
Lrow = zeros(1, Ny*Nx);

for k = 1:size(ray,1)-1
    x1 = ray(k,1);   y1 = ray(k,2);
    x2 = ray(k+1,1); y2 = ray(k+1,2);

    ds = hypot(x2 - x1, y2 - y1);
    if ds == 0, continue; end

    xm = 0.5*(x1 + x2);
    ym = 0.5*(y1 + y2);

    j = floor((xm - x0)/h) + 1;   % x -> column j
    i = floor((ym - y0)/h) + 1;   % y -> row i

    i = max(min(i, Ny), 1);
    j = max(min(j, Nx), 1);

    idx = sub2ind([Ny Nx], i, j);
    Lrow(idx) = Lrow(idx) + ds;
end
end