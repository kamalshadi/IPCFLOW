function [CSb,COy,CO]=cs_ana(CS,sampfreq,maxfreq,segleng,segshift,chanpar);



if chanpar>0
 COy=coherence_partial_complex_b(CS,chanpar);
else
 COy=coherence_plain_complex_b(CS);
end;

CO=abs(COy);

CSb=CS;


return; 