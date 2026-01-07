function [outReal] = makeFFTSymm(inputFFT,fVec,fSamp)
% this function was written to make the fft symmetric and then perform ifft
% to get the real output after gaussian bandpass filter
fMidInd = find(fVec>=fSamp/2,1,'first');
Z1 = inputFFT(1:fMidInd,1);
Z2 = conj(flipud(inputFFT(2:fMidInd,1)));
outReal = ifft([Z1;Z2]);
end

