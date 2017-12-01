function CO=base_coh(CO);
[highcut,nchan,nsegs,nchan]=size(CO);

korr=1+sqrt(eps);
base=CO(:,:,1,:);

%disp([CO(:,47,10,53) ,base(:,47,1,53) ])


for i=1:nsegs
    CO_abs=tanh( atanh(abs(CO(:,:,i,:))/korr)-atanh(abs(base)/korr) );
    CO(:,:,i,:)=CO_abs;
end;


return; 
