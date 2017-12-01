function [ f, val ] = single_sided_fft( X,Fs)
% Fs is sampling frequency   
% X is one dimensional signal
L = length(X); % Length of signal
Y = fft(X);
P2 = Y/L;
val = P2(1:L/2+1);
val(2:end-1) = 2*val(2:end-1);
f = Fs*(0:(L/2))/L;
end
