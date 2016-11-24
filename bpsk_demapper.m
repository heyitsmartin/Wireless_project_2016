function b = bpsk_demapper(a)
% Convert noisy BPSK symbols into a bit vector. Hard decisions.

a = a(:); % Make sure "a" is a column vector

b = real(a) > 0;
b = double(b); % Convert data type from logical to double
end