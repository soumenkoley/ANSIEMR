function [PWSTrace,tVec] = PWS(corrMat,linStack,fSamp,K,nuStack,lagV)
    %This function implemets Phase weighted stacking (PWS) as stated in
    %Schimmel, 2011 (GJI) paper
    
    % Inputs
    % corrMat is a matrix where along each column you have the correlation
    % time series for each time window
    % linStack is the linear stack that will be downweighted
    % fSamp is the sampling frequency in Hz
    % K specifies the width of Gaussian filter, typically Schimmel et al
    % 2011 states K=2
    % nuStack is the power of stacking and its default usage is 2
    
    % Outputs
    % PWSTrace is the phase weighted stacked trace
    
    dt = 1/fSamp;
    [nSamp,nTrace] = size(corrMat);
    tVec = lagV/fSamp;
    
    % first get the S-transform of the linear stack
    [STLin,fVec] = STForward(linStack,dt,K);
    
    % get the phase weights, c_{ps} is Schimmel paper
    sumr = zeros(nSamp-1,nSamp);
    for ii=1:nTrace
        % S-transform of ii-th seismogram
        [stran,fVec] = STForward(corrMat(:,ii),dt,K);
        % Equation 6 in Schimmel et al. (2010; GJI)
%         figure(203)
%         imagesc(tVec,fVec,abs(stran));
%         colormap('jet');
%         colorbar;
        sumr = sumr + (stran./abs(stran)).*exp(1i*2*pi*(fVec'*tVec));
        %disp('I am here');
        
    end
    sumr(:,end) = 0;
    % Equation 6 in Schimmel et al. (2010; GJI)
    sumr = abs(sumr/nTrace).^nuStack;
    
    % apply phase weight to linear stack to get phase-weighted-stack
    % Equation 7 in Schimmel et al. (2010; GJI)
    stranpws = sumr.*STLin;
    
    % transform back to the time domain
    PWSTrace = STInverse(stranpws,fVec);
end

