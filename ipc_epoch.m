function [ f, t, val ] = ipc_epoch(X1,X2,Fs,overlap,seglen,epleng)
% Fs is sampling frequency   
% X is one dimensional signal
L = length(X1);
n_epoch = floor(L/epleng);
'Epoch length'
n_epoch
for j=1:n_epoch
    i1 = (j-1)*epleng+1;
    i2 = i1+epleng-1;
    z1 = X1(i1:i2);
    z2 = X2(i1:i2);
    [ f1, t1, fval1 ] = ipc_t(z1,z2,Fs,overlap,seglen);
    if j==1
        f = f1;
        t = t1;
        n1 = size(fval1,1);
        n2 = size(fval1,2);
        fval = zeros(n1,n2,n_epoch);
        fval(:,:,j) = fval1; 
    else
        fval(:,:,j) = fval1;
    end
end
val = mean(fval,3);
end