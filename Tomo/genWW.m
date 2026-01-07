function [timeWindow] = genWW(vMin,vMax,distOut,fSamp,tWin)
                      % sampling frequency in Hz
    dt = 1/fSamp;                % 0 to 20 s window
    nPts = length(tWin);
    % Compute arrival time bounds
    tMin = distOut/vMax;
    tMax = distOut/vMin;

    % Adjust full time range if needed
    tStart = min(min(tWin), tMin);       % give a bit of buffer if tmin < 0
    tEnd   = max(max(tWin), tMax);      % buffer if tmax > 20

    % New full time vector
    tFull = tStart:dt:tEnd;
    %nFull = length(tFull);

    % Create time-domain window mask
    winMask = (tFull >= tMin) & (tFull <= tMax);  % logical mask
    timeWindow = double(winMask);                 % 1s inside window, 0s outside

    % (Optional) Smooth window edges (taper)
    taper = tukeywin(sum(winMask), 0.2);           % 30% cosine taper
    timeWindow(winMask) = taper;
    timeWindow = timeWindow(1:nPts);
end

