function CS=base_corr(CS);
[highcut,nchan,nsegs,nchan]=size(CS);

base=CS(:,:,1,:);

for i=1:nsegs
    CS(:,:,i,:)=CS(:,:,i,:)./base;
end;

return; 
