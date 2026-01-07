function [S,FVEC] = STForward(s,dt,k)
%
%
% This routine follows from Schimmel, 2011 paper for performing the
% Stockwell transform which is used in Phase weighted stacking
% 
% USAGE: [S,FVEC] = S_transform_FD_fullspec(st,dt,k)
% 
% INPUT:
%   s    = time series
%   dt   = sample interval (s) (default=1)
%   k    = integer value for number of periods to make width of Gaussian
%   (default=2)
% OUPUT:
%   S    = time-frequency matrix of the complex S-transform coefficients
%           (row=f,col=t)
%   FVEC = frequency vector of associated frequencies (symmetric containing
%           both positive and negative frequencies
% 

% AUTHOR:
% Soumen Koley, skoley@nikhef.nl, May 2019
%
% Code mostly from an example by Robert Glenn Stockwell's function
% st.m (beta test)
%--------------------------------------------------------------------------
% set defaults
% making a change to row vector for time series
s = s';
if nargin < 2
    dt = 1; % (s)
    k  = 2; % default is for 2 period window
elseif nargin < 3
    k  = 2; % default is for 2 period window
end
%--------------------------------------------------------------------------
npts=numel(s); % number of points in trace
% make an even number of points
oddflag = 0;
if mod(npts,2)~=0
    s       = s(1:end-1);
    npts    = numel(s);
    oddflag = 1;
end

mid = npts/2;
idx = [0:mid-1, -mid:-1];

FVEC = idx./dt/npts;     % frequency vector
S    = zeros(npts,npts); % allocate S-transform matrix

st = fft(s);     % trace FFT
st = [st st st]; % concatenate trace FFT for sliding to avoid self-aliasing

% loop through frequencies using parallel or serial loops
%
%if matlabpool('size')==0
%     disp('Running S-transform in serial');
    for f=1:npts
        % moving Gaussian window function
        W = exp( -2*pi^2*(idx.^2)./ ( (k/2)*(f-mid-1).^2 ) );
        % slide through FFT and apply window
        S(f,:) = ifft( st(f+mid:f+mid+npts-1).*W);
    end
%else
%     disp('Running S-transform in paralell');
%     parfor f=1:npts
%         % moving Gaussian window function
%         W = exp( -2*pi^2*(idx.^2)./ ( (k/2)*(f-mid-1).^2 ) );
%         % slide through FFT and apply window
%         S(f,:) = ifft( st(f+mid:f+mid+npts-1).*W);
%     end
% end

% this is the f=0 s-transform which is just the mean.
S(mid+1,:) = ifft( st(npts+1:2*npts).*[1 zeros(1,npts-1)] );

if oddflag % add zeros to last time sample is odd npts
    S(:,npts+1) = complex(zeros(npts,1));
end

FVEC=fftshift(FVEC); % arrange frequency vector to be symmetric

%--------------------------------------------------------------------------
return