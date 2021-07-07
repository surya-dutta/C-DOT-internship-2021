function y = convolutional_coding(stream, g)
% INPUT: 
%   bit_stream = uncoded bit stream vector 
%   G = generator matrix 
% OUTPUT: 
%   y = channel coded bit stream 

y = conv2(stream, g);
y = rem(y, 2);
[row, col] = size(y);
y = reshape (y, 1, row * col);
end