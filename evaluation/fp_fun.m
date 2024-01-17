%***********************************************************************
%%% Fletcher-Powell type test function having random local minima and an
%%% arbitrary global minimum set by alfa
clear all

m=5; % Nr. of params.
n=6; %Nr. of extremums 2^n
step=0.5; % resolution


a=50+14.2*randn(n,m);
b=50+14.2*randn(n,m);
alfa=pi/3*randn(m,1);

x=zeros(ceil(2*pi/step),m);
for ijj=1:m
x(:,ijj)=-pi():step:pi();
end


lengx=length(x(:,1));


vec=zeros(lengx,m);
nr=lengx^m;
bb=zeros(n,nr,m);
A=zeros(n,1);
B=zeros(n,nr);

powell_dim=zeros(m,1);
powell_dim(:)=lengx;

fletcher=zeros(powell_dim');
VAR=zeros(nr,m);
IND=zeros(nr,m);
for kk2=1:m
for jj2=1:lengx
IND(jj2*nr/(lengx^kk2)-nr/(lengx^kk2)+1:jj2*nr/(lengx^kk2),kk2)=jj2;
end
 for elem=1:nr
    if IND(elem,kk2)==0
       IND(:,kk2)=repmat(IND(1:elem-1,kk2),lengx^(kk2-1),1);
    end      
 end

end

for ii=1:n
    A(ii)=(a(ii,:)*sin(alfa)+b(ii,:)*cos(alfa));
for kk=1:m
for jj=1:lengx

vec(jj,kk)=a(ii,kk)*sin(x(jj,kk))+b(ii,kk)*cos(x(jj,kk));

end
 for elemc=1:nr
       bb(ii,elemc,kk)=vec(IND(elemc,kk),kk);      
 end

end
  
end

     B(:,:)=sum(bb,3);
     

    for elemf=1:nr
     ind_pointer=num2cell(IND(elemf,:));     
     indx2=sub2ind(size(fletcher),ind_pointer{:});
     %B(ii,elemf)=sum(bb(indx+ii),3);
     fletcher(indx2)=sum(((A(:)-B(:,elemf)).^2)); 
    end
 

    MINF=(min(min(min(min(min(fletcher))))));
    for elemf=1:nr
     ind_pointer=num2cell(IND(elemf,:));     
     indx=sub2ind(size(fletcher),ind_pointer{:}); 
     if fletcher(indx)==MINF
         ind_pointer, elemf1=elemf, x(IND(elemf1,1),1), x(IND(elemf1,2),2), x(IND(elemf1,3),3), x(IND(elemf1,4),4), x(IND(elemf1,5),5)
         alfa
     end
    end

    
% figure(1)
% surf(fletcher)
% hold on
% plot(round((pi+alfa(1))/step),round((pi+alfa(2))/step),'ok')
% hold on
% plot(IND(elemf1,1),IND(elemf1,2),'ok')

minimum=round((pi+alfa)/step)


%save 'fletcher.txt' fletcher -ascii;
% MINF=(min(min(fletcher)));
% 
% for jjj=1:length(x2)
%  for kkk=1:length(x1)
%  if fletcher(jjj,kkk)==MINF
%      jjj,kkk
%      x1(jjj)
%      x2(kkk)
%      alfa
%      clf
%      figure(2)
%      plot(jjj,kkk,'ok')
%      
%  end
%  end
% end
%=======================================================================================================================

% figure(1)
% surf(fletcher)
% hold on 
% plot(round(10*(alfa(1)+pi)+1),round(10*(alfa(2)+pi)+1),'ok')
% 
% MINF
% fletcher(round(10*(alfa(1)+pi)+1),round(10*(alfa(2)+pi)+1))
