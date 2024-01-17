function MS=shiftvec(M,d)

L=max(M);
LL=length(M);
for ii=1:LL
    if M(ii)<=d
       MS(ii)=L-d+M(ii); 
    else
       MS(ii)=M(ii)-d;
    end
end