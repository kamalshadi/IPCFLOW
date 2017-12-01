

function CO=coherence_plain_complex_b(CSal_t)

         
    sqrteps=sqrt(eps);
      
        
        [highcut nchan nsegs nchan]=size(CSal_t);
     
    
     
           
        %Coherence 
          %definition 
            x=reshape(CSal_t,highcut,nchan,nsegs,nchan);
            xout=0*x;

            for j=1:highcut
                for k=1:nsegs
                    y=reshape(x(j,:,k,:),nchan,nchan);
                    z=y./sqrt(diag(y)*(diag(y))'+sqrteps);
                    xout(j,:,k,:)=z;
                end
            end

            clear x y z;
            CO=xout;
            clear xout;

     
       
 return; 


