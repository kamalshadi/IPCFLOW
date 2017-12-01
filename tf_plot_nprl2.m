function [timePoints,freqPoints,finalData] = tf_plot_nprl(data,maxfreq,t0inms,segleng,segshift,titlename)

% makes time-frequency plot   
% usage: tf_plot(data,maxfreq,t0inms,segleng,segshift,titlename)
% 
% By Guido Nolte 
%
finalData = data;
[a,b]=size(data);
clear x y;
for i=1:a
    for j=1:b
        timePoints(i,j)=(i-1)*segshift+segleng/2-t0inms;
        freqPoints(i,j)=(j-1)*maxfreq/(b-1);
    end
end


contourf(timePoints,freqPoints,data,40,'LineColor','none');           % filled contour plot with 40 contour levels, normalized using x and y
shading interp;                       % Each mesh line segment and face has a constant color
maxaxis=max(abs(caxis));caxis([-maxaxis maxaxis]);          % symmetrizes color-axis 
colormap jet;
colorbar                            % Displays colorbar showing the color scale 
title(titlename);   


