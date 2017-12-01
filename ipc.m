function [ f, val ] = ipc(X1,X2,Fs)
% Fs is sampling frequency   
% X is one dimensional signal
[ unused, val1 ] = single_sided_fft( X1,Fs);
[ f, val2 ] = single_sided_fft( X2,Fs);
tmp = val1.*conj(val2);
val = tmp;
end