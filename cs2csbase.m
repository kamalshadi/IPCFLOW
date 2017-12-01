function co_out=cs2csbase(co_in,segs);
% usage co_out=cs2csbase(co_in,segs);

[nfreq,nchan,nseg,nchan]=size(co_in);
[dummy,nsegref]=size(segs);

co_ref=zeros(nfreq,nchan,1,nchan);

for i=1:nsegref
    co_ref=co_ref+co_in(:,:,segs(i),:)/nsegref;
end
%add new line
co_out = zeros(nfreq,nchan,nseg,nchan);
for i=1:nseg
    co_out(:,:,i,:)=co_in(:,:,i,:)-co_ref;
end

return;



