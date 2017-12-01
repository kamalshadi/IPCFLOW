function [CSbase,COy,CO]=cs_ana(CSal,sampfreq,maxfreq,segleng,segshift,chanpar);
%renamed CSb to CSbase and CS to CSal to match main function


if chanpar>0
 COy=coherence_partial_complex_b(CSal,chanpar);
else
 COy=coherence_plain_complex_b(CSal);
end;

CO=base_coh(COy);

COy=base_coy(COy);
CSbase=base_corr(CSal);


return; 