function B=mirror(A,side,mode,domdef)
domx=length(A(:,1));
domy=length(A(1,:));
B=zeros(size(A));
if strcmp(side,'row')==1
for ii=1:domx
    B(ii,:)=A(domx+1-ii,:);
end
elseif strcmp(side,'column')==1
for ii=1:domy
    B(:,ii)=A(:,domy+1-ii);
end
else
for ii=1:domx
    B(ii,:)=A(domx+1-ii,:);
end
for ii=1:domy
    B(:,ii)=A(:,domy+1-ii);
end   
end

if mode=='coord' 
    domx=domdef;
    if strcmp(side,'row')==1
    B(1,:)=A(1,:);
    B(2,:)=domx+1-A(2,:);
    B(3,:)=A(3,:);
    end
end
end