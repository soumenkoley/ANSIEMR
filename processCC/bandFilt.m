function [filtSig] = bandFilt(inpSig,fSamp,alpha,beta,fC)
    %This function applies a narrow pass band filter to the data
    lenSig = length(inpSig);
    fftSig = fft(inpSig,lenSig);
    if(rem(lenSig,2)==0)
        fftSigUse = fftSig(1:1:(lenSig/2+1),1);
        fSig = fSamp/2*linspace(0,1,(lenSig/2+1));
    else
        fftSigUse = fftSig(1:1:((lenSig-1)/2+1),1);
        fSig = fSamp/2*linspace(0,1,((lenSig-1)/2+1));
    end
    
    [filterVal] = getGaussFilter(fC,alpha,fSig,beta);
    filtFFT = fftSigUse.*filterVal';
    
    symFFT = make_FFT_symmetricNew(filtFFT);
    
    filtSig = ifft(symFFT);
    filtSig = [filtSig;filtSig(end,1)];
end