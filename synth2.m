function z = synth2(fs,duration,f1,f2,s,da,mod)
% this function synthesize a modulated signal (f1 modulate f2)
% fs is sampling rate
% T is dutation
% f1 is the lower frequecy
% f2 is higher frequency
% s is the snr
ts = 1/fs;
t = 0:ts:(duration-ts);
delay = 2*pi*f2*da*0.001;
modulated_signal = t;
signal = t;
for i=1:length(t)
    signal(i) = sqrt(2)*cos(2*pi*f1*t(i));
    if (mod)
        modulated_signal(i) = sqrt(2)*cos(2*pi*f2*t(i)-delay).*signal(i);
    else
        modulated_signal(i) = sqrt(2)*cos(2*pi*f2*t(i)-delay);
    end
end
z1 = awgn(modulated_signal,s);
z2 = awgn(signal,s);
z = [z1;z2]';
% plot(t(1:5*fs),z(1:5*fs));
% N = length(z);
% xdft = fft(z);
% xdft = xdft(1:N/2+1);
% psdx = (1/(fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:fs/N:fs/2;
% figure;
% plot(freq,psdx)
% figure
% spectrogram(modulated_signal,128,5,[1:40],fs,'yaxis')
end