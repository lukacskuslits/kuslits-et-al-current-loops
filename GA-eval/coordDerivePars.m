function coord=coordDerivePars(pd_t,h,thres)
s=size(pd_t);
kk=1;

def=min(min(pd_t));
for ii=2:s(1)-1
    for jj=2:s(2)-1
        if (pd_t(ii,jj)>def+h && pd_t(ii-1,jj)>def+h)...
           || (pd_t(ii,jj)>def+h && pd_t(ii+1,jj)>def+h)...
           || (pd_t(ii,jj)>def+h && pd_t(ii,jj-1)>def+h)...
           || (pd_t(ii,jj)>def+h && pd_t(ii,jj+1)>def+h)
            if abs(gradient(pd_t(ii,jj)))<thres...
               && pd_t(ii,jj)==max(max(pd_t(ii-1:ii+1,jj-1:jj+1)))
                coord(:,kk)=[ii;jj];
                kk=kk+1;
            end
        end
    end
end
% zed=10*ones(1,length(coord(1,:)));
% figure(fig), hold on, scatter3(coord(2,:),coord(1,:),zed,300,'black','filled')
% ylim([1 45])
end