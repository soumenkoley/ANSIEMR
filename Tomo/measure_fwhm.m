function fwhm_t = measure_fwhm(col, tvec, pk_idx)
% measure_fwhm - measure full width at half maximum (FWHM) of FTAN peak
%
% Inputs:
%   col    - FTAN amplitude column (vector over time)
%   tvec   - corresponding time axis (same length as col)
%   pk_idx - index of picked maximum in col
%
% Output:
%   fwhm_t - FWHM in seconds (time units of tvec).
%            Returns NaN if undefined (flat or zero peak).

    % normalize column to its maximum
    amp = col(:) / max(col);
    peak_val = amp(pk_idx);
    if peak_val <= 0
        fwhm_t = NaN;
        return;
    end

    half_height = 0.5 * peak_val;

    % search left
    i_left = pk_idx;
    while i_left > 1 && amp(i_left) > half_height
        i_left = i_left - 1;
    end
    if i_left == 1
        t_left = tvec(1);
    else
        % linear interpolate between i_left and i_left+1
        frac = (half_height - amp(i_left)) / (amp(i_left+1) - amp(i_left) + eps);
        t_left = tvec(i_left) + frac*(tvec(i_left+1)-tvec(i_left));
    end

    % search right
    i_right = pk_idx;
    while i_right < numel(amp) && amp(i_right) > half_height
        i_right = i_right + 1;
    end
    if i_right == numel(amp)
        t_right = tvec(end);
    else
        frac = (half_height - amp(i_right)) / (amp(i_right-1) - amp(i_right) + eps);
        t_right = tvec(i_right) + frac*(tvec(i_right-1)-tvec(i_right));
    end

    fwhm_t = abs(t_right - t_left);
end