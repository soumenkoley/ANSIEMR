function N = make_FFT_symmetricNew(Z)
% Make Z symmetric in the sense that the FFT of the output N will be real

% S. Koley, Nikhef
% April, 2016

[nRows, nCols] = size(Z);

%evenNrows = mod(nRows, 2) == 0;

N = zeros((nRows-1)*2,nCols);

N(1:nRows, :) = Z(1:nRows, :);
N(nRows+1:((nRows-1)*2), :) = conj( flipud1( N(2:(nRows-1), :) ) );

end

function y = flipud1(x)

y = x(end:-1:1,:); 

end

