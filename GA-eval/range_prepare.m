function [range]=range_prepare(l_array,h_array)
rl_array=l_array';
rh_array=h_array';
range=zeros(2*length(rl_array(:,1)),length(rl_array(1,:)));
jj=1; kk=1;
for ii=1:length(range(:,1))
    if mod(ii,2)==1
        range(ii,:)=rl_array(jj,:);
        jj=jj+1;
    end
    if mod(ii,2)==0
        range(ii,:)=rh_array(kk,:);
        kk=kk+1;
    end
end

for ll=1:length(range(1,:))
for ii=1:2:length(range(:,1))
if range(ii+1,ll)<range(ii,ll)
    range(ii:ii+1,ll)=swaprows(range(ii:ii+1,ll),1,2);
end
end
end
end
    