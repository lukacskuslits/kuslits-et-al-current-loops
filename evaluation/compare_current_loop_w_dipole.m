clear all, close all
mu = 4*pi*1e-7;
lat = 0;
lon = 0;
Roc = 3.48e6;
z = 1e5;
depth = Roc - z;
r = 1e6;
I = 1e9;
t = lat;
l = lon;
attenuate = 1100;
rmin = 5e5;
depth_max = 3.48e6; %3.48e6;
depth_min = 2.7e6;
dBt_max = 7.3e-13; %4e-11 with the lowest std
pp0 = [lat;lon;r;I;t;l;depth;0];
pp0(5,:)=pp0(1,:)*pi()/180;
pp0(6,:)=pp0(2,:)*pi()/180;
pp0(1,:)=pp0(1,:)*pi()/180;
pp0(2,:)=pp0(2,:)*pi()/180;
res_loop = total_field_approx(pp0, 100, 200);
res_loop = reshape(res_loop,200,100)';

% res_dip = total_field_approx_dip(pp0,100,200);
% res_dip = reshape(res_dip,200,100)';

% plot(res_dip(:,100)*4*pi,'r')
% hold on
plot(res_loop(:,100))

save('res_dip_1e6.mat','res_dip')
save('res_loop_1e6.mat','res_loop')
