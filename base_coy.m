function CO=base_coy(CO);
[highcut,nchan,nsegs,nchan]=size(CO);

korr=1+sqrt(eps);
base=CO(:,:,1,:);

%disp([CO(:,47,10,53) ,base(:,47,1,53) ])


for i=1:nsegs
    CO_r=tanh( atanh(real(CO(:,:,i,:))/korr)-atanh(real(base)/korr) );
    CO_i=tanh( atanh(imag(CO(:,:,i,:))/korr)-atanh(imag(base)/korr) );
    CO(:,:,i,:)=CO_r+sqrt(-1)*CO_i;
   if i==1
       %disp([CO(1,1,1,2) ,base(1,1,1,2) ])
   end
end;

return; 
