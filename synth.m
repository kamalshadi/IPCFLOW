function z = synth(fs,ncycle,f1,f2,s)
% this function synthesize a modulated signal (f1 modulate f2)
% fs is sampling rate
% T is dutation
% f1 is the lower frequecy
% f2 is higher frequency
% s is the snr
ts = 1/fs;
t = 0:ts:(ncycle-1)*fs;
high_frequency_signal = sqrt(2)*cos(2*pi*f2*t);
low_frequency_signal = sqrt(2)*cos(2*pi*f1*t);
modulated_signal = high_frequency_signal.*low_frequency_signal+1;
z = awgn(modulated_signal,s);
plot(t(1:5*fs),z(1:5*fs));
N = length(z);
xdft = fft(z);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/N:fs/2;
figure;
plot(freq,psdx)
figure
spectrogram(modulated_signal,128,5,[1:40],fs,'yaxis')
end