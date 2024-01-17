function Fmvj=Fmpart(pnew,subname,m,looppar,n,scale,domx,domy)
PN=zeros(m,n);
 for i=1:n
       PN(1:m,i)=pnew(m*i-m+1:m*i);
       %pnew(:,i)=pN(looppar*i-looppar+1:looppar*i);
    Fmvtemp=log(abs(subname(PN,looppar,n,scale,domx,domy)));
 end
 for i=1:domy
 Fmvj(domx*i-domx+1:domx*i)=Fmvtemp(:,i);
 end
end