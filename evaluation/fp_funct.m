%***********************************************************************
%%% Fletcher-Powell type test function having random local minima and an
%%% arbitrary global minimum set by alfa
function [target]=fp_funct(pf,looppar,f,scale,alfa)

%looppar=4; % Nr. of params.
n=6; %Nr. of extremums 2^n
%step=0.5; % resolution
%=100;

a=50+14.2*randn(n,looppar);
b=50+14.2*randn(n,looppar);


x=pf(:,1,1,1,1)/scale;


% x=zeros(looppar,nloop,1);
% for nloops=1:nloop
% 
% x(:,nloops)=pf(:,nloops,1,1,1)'/scale;
% 
% end


bb=zeros(n,looppar);
A=zeros(n,1);


for ii=1:n
    A(ii)=(a(ii,:)*sin(alfa)+b(ii,:)*cos(alfa));
for kk=1:looppar

bb(ii,kk)=a(ii,kk)*sin(x(kk))+b(ii,kk)*cos(x(kk));

end    

end
  
B=sum(bb,2);


fletcher=sum(((A(:)-B(:)).^2));     

target=1/abs(fletcher-f);

 
       end
       
%     MINF=(min(min(min(min(min(fletcher))))));
%     for elemf=1:nr
%      ind_pointer=num2cell(IND(elemf,:));     
%      indx=sub2ind(size(fletcher),ind_pointer{:}); 
%      if fletcher(indx)==MINF
%          ind_pointer, elemf1=elemf, x(IND(elemf1,1),1), x(IND(elemf1,2),2), x(IND(elemf1,3),3), x(IND(elemf1,4),4), x(IND(elemf1,5),5)
%          alfa
%      end
%     end

    
% figure(1)
% surf(fletcher)
% hold on
% plot(round((pi+alfa(1))/step),round((pi+alfa(2))/step),'ok')
% hold on
% plot(IND(elemf1,1),IND(elemf1,2),'ok')




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
