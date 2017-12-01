%function plot_coh_onepair(CO,sampfreq,segleng,segshift,t0,chan1,chan2,chan_names,xtitle);
function [timePoints,freqPoints,finalData] = plot_coh_onepair(CO,sampfreq,segleng,segshift,t0,chan1,chan2);

[highcut,nchan,nsegs,nchan]=size(CO);

T=segleng/sampfreq;
maxfreqexact=(highcut-1)/T;
t0inms=t0*1000/sampfreq;

%[n,m]=size(locs);
%ind2chan=locs(:,1);
%nplot=nc;

%figure;
   %for i=1:1;for j=1:1;
       CO_here=(reshape(CO(:,chan1,:,chan2),highcut,nsegs))';
       %titlename=strcat(xtitle,chan_names(chan1),',',chan_names(chan2));
       %subplot(nplot-1,nplot-1,(i-2)*(nplot-1)+j);
       [timePoints,freqPoints,finalData] = tf_plot_nprl2(CO_here,maxfreqexact,t0inms,segleng*1000/sampfreq,segshift*1000/sampfreq,'');
       xlabel('Time [ms]');
       ylabel('Frequency [Hz]');
       %  end;end;


return;