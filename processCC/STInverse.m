function st = STInverse(S,fvec)
%
% This script calculates the inverse S-transfrom of a signal by integrating
% over the time window and inverse Fourier transforming.
%
% USAGE: st = S_transform_inverse_fullspec(S,fvec)
%
% INPUT:
%   S    = time-frequency matrix of the complex S-transform coefficients
%           (row=f,col=t)
%   fvec = Frequency vector (symmetric contain both positive and negative
%           frequencies)
% OUPUT:
%   st   = time series trace after inverse S-transform
%
%
% AUTHOR:
% Soumen Koley, skoley@nikhef.nl, May 2019
%
% Code mostly from an example by Robert Glenn Stockwell's function
% inverse_st.m ($Revision: 1.0 $  $Date: 2004/10/10  $)
%
% Reference is "Localization of the Complex Spectrum: The S-Transform" from
% IEEE Transactions on Signal Processing, vol. 44., number 4, April 1996,
% pages 998-1001.

% Dimensions of S-transform matrix
[nf,npts] = size(S);

% Check that matrix is in correct dimensions
if nf ~= numel(fvec)
    S = S';
    [~,npts] = size(S);
end

% integrate over time-axis to get FFT spectrum
spec_full = sum(S,2); 

% the time series is the inverse fft of this
ts = ifft(fftshift(spec_full));

% and take the real part, the imaginary part should be zero.
st = real(ts);

% if odd number of points add 1 sample to match original trace
if rem(npts,2) ~= 0
    st(end+1)=0;
end

return