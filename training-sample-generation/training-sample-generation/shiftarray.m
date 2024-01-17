function MS=shiftarray(M,d,side)


%disp(M)

% % if strcmp(mode,'grid')==1 
% %     if strcmp(side,'row')==1
%     switch side
%         case 'row'
%         case 'col'
%         MS(2,:)=shiftvec(M(2,:),d);
%         MS(1,:)=M(1,:);
%     end
%       %B(3,:)=M(3,:);
% %     end
% % else
        SM=size(M);
    switch side
    
    case 'row'
    case 'col'
       V=1:SM(2);
       VS=shiftvec(V,d);
       
       MS(:,V)=M(:,VS);
    end
end
    
   