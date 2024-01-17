function [l_array,h_array]=lh_array_new(array, abs, domx, domy)
l_array=array-abs;
h_array=array+abs;
h_array(h_array>domy) = domy;
hlat = h_array(1:2:end,:);
hlat(hlat>domx) = domx;
h_array(1:2:end,:) = hlat;
l_array(l_array<1)=1;
end