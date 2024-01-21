function MS=shiftvec(M,d)
%% Performs translation on a vector
%Input parameters:
%-----------------
% M: input vector
% d: translation distance (elements)
%Ouptut values:
%--------------
% MS: shifted (translated) vector

L=max(M);
LL=length(M);
for ii=1:LL
    if M(ii)<=d
       MS(ii)=L-d+M(ii); 
    else
       MS(ii)=M(ii)-d;
    end
end
