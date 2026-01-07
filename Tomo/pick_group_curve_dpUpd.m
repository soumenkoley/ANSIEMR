function out = pick_group_curve_dpUpd(ftan, freqs, times, distance_m, varargin)
% PICK_GROUP_CURVE_DP  Pick smooth group-velocity curve from FTAN via DP
%
% Inputs:
%   ftan       : matrix (n_freq x n_time)
%   freqs      : vector length n_freq (Hz)
%   times      : vector length n_time (s)
%   distance_m : scalar, inter-station distance (m)
%
% Optional name-value inputs:
%   'vel_axis', 'vmin', 'vmax', 'smoothness_lambda', 'max_jump_m_s',
%   'weight_peak_boost', 'peak_prominence', 'energy_exponent',
%   'amp_floor', 'post_sg_win', 'post_sg_poly', 'debug'
%
% Output struct fields:
%   times_pick, vel_pick, slowness_pick, cost_map, path_cost, mask_valid, idx_path

% ---- parse inputs ------------------------------------------------------
p = inputParser;
addRequired(p, 'ftan');
addRequired(p, 'freqs');
addRequired(p, 'times');
addRequired(p, 'distance_m');
addParameter(p, 'vel_axis', [], @(x) isempty(x) || isvector(x));
addParameter(p, 'vmin', 300, @isnumeric);
addParameter(p, 'vmax', 5000, @isnumeric);
addParameter(p, 'smoothness_lambda', 1e6, @isnumeric);
addParameter(p, 'max_jump_m_s', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'weight_peak_boost', true, @islogical);
addParameter(p, 'peak_prominence', 0.05, @isnumeric);
addParameter(p, 'energy_exponent', 1.0, @isnumeric);
addParameter(p, 'amp_floor', 1e-12, @isnumeric);
addParameter(p, 'post_sg_win', 11, @isnumeric);
addParameter(p, 'post_sg_poly', 2, @isnumeric);
addParameter(p, 'debug', false, @islogical);
parse(p, ftan, freqs, times, distance_m, varargin{:});
opt = p.Results;

[n_freq, n_time] = size(ftan);
if length(freqs) ~= n_freq, error('freqs must match ftan rows'); end
if length(times) ~= n_time, error('times must match ftan cols'); end

% ---- velocity axis ------------------------------------------------------
if isempty(opt.vel_axis)
    vel_axis = nan(1, n_time);
    pos_idx = times > 0;
    vel_axis(pos_idx) = distance_m ./ times(pos_idx);
    vel_axis(~pos_idx) = Inf;
else
    vel_axis = opt.vel_axis(:).'; % row vector
end
slowness = 1 ./ vel_axis;
slowness(~isfinite(slowness)) = Inf;

% ---- preprocessing & normalization -------------------------------------
ftan = double(ftan);
ftan(ftan < 0) = 0;

% Gaussian smoothing along time axis
ftan_sm = ftan;
for i = 1:n_freq
    ftan_sm(i,:) = smoothdata(ftan(i,:), 'gaussian', 5);
end

% normalize rows
row_max = max(ftan_sm, [], 2);
row_max(row_max < opt.amp_floor) = opt.amp_floor;
ftan_norm = ftan_sm ./ row_max;

% apply energy exponent
if opt.energy_exponent ~= 1.0
    ftan_norm = ftan_norm .^ opt.energy_exponent;
end

% optional peak boosting
peak_boost = ones(size(ftan_norm));
if opt.weight_peak_boost
    for i = 1:n_freq
        col = ftan_norm(i,:);
        prom = opt.peak_prominence * (max(col)+eps);
        [~, locs] = findpeaks(col, 'MinPeakProminence', prom);
        if ~isempty(locs)
            peak_boost(i, locs) = 2;
        end
    end
    ftan_norm = ftan_norm .* peak_boost;
end

% ---- velocity mask ------------------------------------------------------
if numel(opt.vmin)==1, vmin_vec = opt.vmin*ones(n_freq,1); else vmin_vec = opt.vmin(:); end
if numel(opt.vmax)==1, vmax_vec = opt.vmax*ones(n_freq,1); else vmax_vec = opt.vmax(:); end

mask_valid = false(n_freq,n_time);
for i = 1:n_freq
    mask_valid(i,:) = (vel_axis >= vmin_vec(i)) & (vel_axis <= vmax_vec(i));
end

% ensure at least one valid starting point
if all(~mask_valid(1,:))
    [~, idx] = min(abs(vel_axis - mean([vmin_vec(1) vmax_vec(1)])));
    mask_valid(1, idx) = true;
end

% --- EDGE REPEL weighting -----------------------------------------------
K = 100; alpha = 0.5; psh = 2;
edge_weight = ones(n_freq, n_time);
for i = 1:n_freq
    vmin = vmin_vec(min(i, numel(vmin_vec)));
    vmax = vmax_vec(min(i, numel(vmax_vec)));
    dv_low  = max(0, vel_axis - vmin);
    dv_high = max(0, vmax - vel_axis);
    d_edge  = min(dv_low, dv_high);
    inside = (vel_axis >= vmin) & (vel_axis <= vmax);
    near   = inside & (d_edge <= K);
    w = 1 - alpha*(1 - (d_edge./K)).^psh;
    edge_weight(i, near) = w(near);
end
ftan_norm = ftan_norm .* edge_weight;

% ---- local cost --------------------------------------------------------
big = 1e300;
local_cost = -ftan_norm;
local_cost(~mask_valid) = big;

% ---- DP arrays ---------------------------------------------------------
cum_cost = big * ones(n_freq, n_time);
backptr  = -ones(n_freq, n_time, 'int32');
cum_cost(1, mask_valid(1,:)) = local_cost(1, mask_valid(1,:));

% ---- DP main loop ------------------------------------------------------
for i = 2:n_freq
    prev_costs = cum_cost(i-1,:);
    valid_prev = isfinite(prev_costs) & (prev_costs < big/2);
    if ~any(valid_prev)
        cum_cost(i,:) = local_cost(i,:); 
        continue;
    end
    
    for r = 1:n_time
        if ~mask_valid(i,r), continue; end
        sl_r = slowness(r);
        sq = (sl_r - slowness).^2;
        trans = opt.smoothness_lambda * sq;
        if ~isempty(opt.max_jump_m_s)
            v_r = vel_axis(r);
            trans(abs(vel_axis - v_r) > opt.max_jump_m_s) = big;
        end
        costs = prev_costs;
        costs(~valid_prev) = big;
        total_costs = costs + trans;
        [best_cost, best_idx] = min(total_costs);
        if isfinite(best_cost) && best_cost < big/2
            cum_cost(i,r) = local_cost(i,r) + best_cost;
            backptr(i,r) = int32(best_idx);
        end
    end
end

% ---- backtracking ------------------------------------------------------
last_col = cum_cost(n_freq,:);
valid_end = isfinite(last_col) & (last_col < big/2);

path = -ones(1,n_freq,'int32');
times_pick = nan(1,n_freq);
vel_pick   = nan(1,n_freq);
sl_pick    = nan(1,n_freq);
fwhm_pick  = nan(1,n_freq);   % new: FWHM (s)
sigma_pick = nan(1,n_freq);   % new: Gaussian sigma equivalent

if any(valid_end)
    [~, end_idx] = min(last_col);
    path(n_freq) = int32(end_idx);
    for i = n_freq:-1:2
        prev = backptr(i, path(i));
        if prev < 1, break; end
        path(i-1) = prev;
    end
    for i = 1:n_freq
        idx = path(i);
        if idx >= 1
            times_pick(i) = times(idx);
            vel_pick(i)   = vel_axis(idx);
            sl_pick(i)    = slowness(idx);
            % --- measure FWHM for uncertainty ---
            col = ftan(i,:);       % FTAN row at this frequency
            tvec = times;          % same axis
            pk_idx = idx;          % index of chosen travel-time
            fwhm = measure_fwhm(col, tvec, pk_idx);
            fwhm_pick(i) = fwhm;
            sigma_pick(i) = fwhm / (2*sqrt(2*log(2)));  % Gaussian ?
        end
    end
end

% ---- post smoothing ---------------------------------------------------
valid = isfinite(sl_pick);
if sum(valid)>=5
    x = 1:n_freq;
    interp_sl = interp1(x(valid), sl_pick(valid), x, 'linear', 'extrap');
    sgw = double(opt.post_sg_win);
    if mod(sgw,2)==0, sgw=sgw+1; end
    if sgw > n_freq, sgw = n_freq-(1-mod(n_freq,2)); end
    try
        sl_sm = sgolayfilt(interp_sl, opt.post_sg_poly, sgw);
        sl_sm(~valid)=NaN;
        vel_sm = 1./sl_sm;
        vel_pick = vel_sm;
        sl_pick = sl_sm;
        if ~isempty(distance_m)
            times_pick = distance_m ./ vel_pick;
        end
    catch
    end
end

% ---- outputs -----------------------------------------------------------
out = struct();
out.times_pick   = times_pick;
out.vel_pick     = vel_pick;
out.slowness_pick= sl_pick;
out.cost_map     = cum_cost;
out.path_cost    = (path(end) > 0) * cum_cost(n_freq, max(path(end),1));
out.mask_valid   = mask_valid;
out.idx_path     = path;
% new outputs
out.fwhm_pick  = fwhm_pick;
out.sigma_pick = sigma_pick;

if opt.debug
    out.ftan_norm = ftan_norm;
    out.peak_boost = peak_boost;
end
end