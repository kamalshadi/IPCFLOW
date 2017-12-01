function [ f_crop, t_crop, X_crop ] = ipc_t_crop(t,f,X,tmin,tmax,fmin,fmax)
i1 = length(t);
i2 = length(f);
tmax_set = 0;
tmin_set = 0;
fmax_set = 0;
fmin_set = 0;
i_tmin = 1;
i_tmax = i1;
i_fmin = 1;
i_fmax = i2;
for i=1:i1
    if (t(i)>=tmin && ~tmin_set) 
        tmin_set=1;i_tmin=i;
    end;
    if (t(i)>tmax && ~tmax_set) 
        tmax_set=1;i_tmax=i-1;
    end;
    if (tmin_set && tmax_set)
        break;
    end;
end
for i=1:i2
    if (f(i)>=fmin && ~fmin_set) 
        fmin_set=1;i_fmin=i;
    end;
    if (f(i)>fmax && ~fmax_set) 
        fmax_set=1;i_fmax=i-1;
    end;
    if (fmin_set && fmax_set)
        break;
    end;
end
f_crop = f(i_fmin:i_fmax);
t_crop = t(i_tmin:i_tmax);
X_crop = X(i_fmin:i_fmax,i_tmin:i_tmax);
end