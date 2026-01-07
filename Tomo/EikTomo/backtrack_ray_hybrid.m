function ray = backtrack_ray_hybrid(T, x0, y0, h, xs, ys, xr, yr, ...
                                   step_frac, max_steps, tol_dist)

[Ny, Nx] = size(T);
ray = zeros(max_steps+2, 2);

x = xr; 
y = yr;
ray(1,:) = [x, y];
prev_T = interp_bilinear(T, x, y, x0, y0, h);

step_len = step_frac * h;
min_step_len = 0.1 * h;

for k = 1:max_steps
    
    %% Continuous gradient descent step
    [Tx, Ty] = gradT_at_point(T, x, y, x0, y0, h);
    gnorm = hypot(Tx, Ty);
    
    if gnorm > 1e-12
        dx = -step_len * Tx / gnorm;
        dy = -step_len * Ty / gnorm;

        x_new = min(max(x + dx, x0), x0 + (Nx-1)*h);
        y_new = min(max(y + dy, y0), y0 + (Ny-1)*h);

        Tnew = interp_bilinear(T, x_new,y_new, x0,y0,h);

        if Tnew < prev_T
            x = x_new; y = y_new;
            prev_T = Tnew;
            ray(k+1,:) = [x,y];

            if hypot(x-xs,y-ys) < tol_dist
                ray(k+2,:) = [xs,ys];
                ray = ray(1:k+2,:);
                return;
            end

            continue; % good step, go to next iteration
        end
    end
    
    %% Fallback: move to best neighboring node
    [i,j] = xy_to_ij(x,y,x0,y0,h,Ny,Nx);
    
    best_i = i; best_j = j; best_T = prev_T;
    for di=-1:1
        for dj=-1:1
            if di==0 && dj==0, continue; end
            ii = i+di; jj = j+dj;
            if ii<1 || ii>Ny || jj<1 || jj>Nx, continue; end
            if T(ii,jj) < best_T
                best_T = T(ii,jj);
                best_i = ii; best_j = jj;
            end
        end
    end
    
    if best_T < prev_T
        x = x0 + (best_j-1)*h;
        y = y0 + (best_i-1)*h;
        prev_T = best_T;
        ray(k+1,:) = [x,y];
        continue;
    end
    
    %% Stop if both methods stuck
    if step_len < min_step_len
        ray(k+1,:) = [xs,ys];
        ray = ray(1:k+1,:);
        return;
    end
    
    step_len = step_len * 0.5; % try smaller step
end

ray(end,:) = [xs,ys];
end