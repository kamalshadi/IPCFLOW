function [ f, t, fval ] = ipc_t(X1,X2,Fs,overlap,seglen)
% Fs is sampling frequency   
% X is one dimensional signal
L = length(X1);
fwd = floor((1-overlap)*seglen);
time_resolution_in_sample =fwd;
num_time_bins = floor((L-seglen)/fwd)+1;
time_resolution_in_sec = fwd/Fs;
Freq_resolution_in_hz = Fs/seglen;
t = zeros(1,num_time_bins);
num_time_bins*time_resolution_in_sec;
t = linspace(0,num_time_bins*time_resolution_in_sec,num_time_bins);
f = 0;
for i=1:num_time_bins
    t(i) = ((i-1)*fwd+seglen/2)/Fs;
    if (i>1)
        i1 = i1+fwd;
        i2 = i1+seglen-1;
    else
        i1 = 1;
        i2 = i1+seglen-1;
    end

%     [num2str(i2) 'from' num2str(L)]
    z1 = X1(i1:i2);
    z2 = X2(i1:i2);
%     Y1 = prctile(z1,90);
%     Y2 = prctile(z2,90);
%     z1(z1>Y1) = Y1;
%     z2(z2>Y2) = Y2;
%     z1(z1<-Y1) = -Y1;
%     z2(z2<-Y2) = -Y2;
    sig1 = z1 - mean(z1);
    sig2 = z2 - mean(z2);
    [tmp,val]=ipc(sig1,sig2,Fs);

    if (i==1) 
        f=tmp; 
        l = length(val); 
        fval = zeros(l,num_time_bins); 
        fval(:,i)=val;
    else
        fval(:,i)=val;
    end
end
end