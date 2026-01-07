function [f, c, mask] = phaseVel(ccf_t, fs, delta_km, fmin, fmax, win, poly)
    n = numel(ccf_t);
    %w = tukeywin(n,0.1);
    C = fft(ccf_t(:));
    freqs = (0:floor(n/2))' * fs / n;
    C = C(1:numel(freqs));
    band = freqs>=fmin & freqs<=fmax;
    f = freqs(band);
    Cb = C(band);

    % Savitzky-Golay smooth real & imag
    %Re = sgolayfilt(real(Cb), poly, win);
    %Im = sgolayfilt(imag(Cb), poly, win);
    %Csm = Re + 1i*Im;
    
    Csm = Cb;
    omega = 2*pi*f;
    dC = gradient(Csm)./gradient(omega); % dC/d?

    amp = abs(Csm);
    eps = 0.05*median(amp);
    invC = conj(Csm) ./ max(amp.^2, eps.^2);

    phi_prime = imag(invC .* dC);
    tau_p = phi_prime;

    c = (delta_km) ./ tau_p;

    cmin = 0.4e3; cmax = 5.0e3;
    mask = amp > prctile(amp,30) & isfinite(c) & c>cmin & c<cmax;

    % smooth slowness on valid points
    sl = nan(size(c)); sl(mask) = 1./c(mask);
    sl = fillmissing(sl,'nearest');
    sl = sgolayfilt(sl, 2, win);
    c = 1./sl;
end