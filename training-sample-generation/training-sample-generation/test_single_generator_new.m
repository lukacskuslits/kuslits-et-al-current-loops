clear all, close all
dt = 25*3.15e7;
mu = 4*pi*1e-7;

% coef = load('coef_new.mat');
% coef = coef.coef;
% powers = load('powers_new.mat');
% powers = powers.powers;
% X = ones(16200,3);
% powers = cast(powers, 'like', coef);
% XX = repmat(X,[1,1,length(powers)]);
% powers = repmat(powers,[1,1, length(XX)]);
% powers = permute(powers, [3,2,1]);
% save('coef_use.mat', 'coef')
% save('powers_use.mat', 'powers')

n = 1;  
rmin = 3.25; %3.25; %original upper and lower vals for R
rmax = 10; %10; %
Imin = 2.2e8;
Imax = 1e9;
depth_min = 2.6; %2.6; ICB
depth_max = 3.4; %3.48; CMB
dtB_max = 7e-13; %4e-11 with the lowest std
%Highest possible rate of change in current:
%dI/dt = 2*Rmax*(dB/dt)max/mu
dtI_max = 2*(((depth_max-depth_min)*1e6)^2+(3.25*1e5)^2)^(3/2)*(mu*(3.25*1e5)^2)^(-1)*dtB_max;
dtI_max = dtI_max/1100;
disp('dIt')
disp(dtI_max)
dtI_min = 0;
v_max = 0;

%% Primary Field
limits_static = [1;89;1;179;rmin;rmax;Imin;Imax;-1;1;0;0;0;0;depth_min;depth_max;dtI_min;0;-v_max;v_max];
pp_static = loop_gen_constr_rad_total(limits_static,n,3,3,3,3,[2,2]); %read_results
disp(pp_static)
%[coef, pow] = transform_coefficients();
%tic
deg_res = 2;
cords=pp_static(1:2,:)/deg_res;
pp_static(5,:)=pp_static(1,:)*pi()/180;
pp_static(6,:)=pp_static(2,:)*pi()/180;
pp_static(1,:)=pp_static(1,:)*pi()/180;
pp_static(2,:)=pp_static(2,:)*pi()/180;
res_static = total_field_approx_new(pp_static, 90, 180);
%toc
stat = reshape(res_static, 180, 90)';
figure(2), imagesc(stat)

%Test positional marks
lat_dim = 180/deg_res;
long_dim = 360/deg_res;
[pmarks,radmarks,Imarks_neg,Imarks_pos,depthmarks,dtImarks_neg,dtImarks_pos,vlats_neg,vlats_pos,vlongs_neg,vlongs_pos]=pos_marks([cords;pp_static(3:end,:)],lat_dim,long_dim,1,1);
v = [5, 5];
h = [5, 5];
pos_data_i = gaussing(pmarks, 3, v, h);
figure(3), imagesc(pos_data_i)
%% Total Field
pp_total = pp_static;
pp_total(8,:) = pp_total(4,:)/2e11;
% res_total = total_field_approx_orig(pp_total, 90, 180);
% tot = reshape(res_total, 180, 90)';
% figure(2), imagesc(tot)
% disp(max(max(abs(tot))))
tic
res_total = total_field_approx_new(pp_total, 90, 180);
toc
tot = reshape(res_total, 180, 90)';
% figure(3), imagesc(tot)
% disp(max(max(abs(tot))))

%TODO: ide tot helyett majd az eredeti ter kellene
%TODO: Hogyan normáljam 0 kozeli helyeken?
ind_field = tot-stat;
% figure(4), imagesc(ind_field)

save('ind_field.mat', 'ind_field')
save('tot_single.mat','tot')
save('stat_single.mat','stat')

%% Plot fields
clear all, close all
load('stat_single.mat')
figure(3), imagesc([0:360],[0:180],stat*1000)
ax = get(gca,'XTick');
set(gca,'XTick',[0:60:360])
set(gca,'GridLineStyle','-')
set(gca,'LineWidth',2)
grid on
colorbar();

load('tot_single.mat')
figure(4), imagesc([0:360],[0:180],tot*1000)
ax = get(gca,'XTick');
set(gca,'XTick',[0:60:360])
set(gca,'GridLineStyle','-')
set(gca,'LineWidth',2)
grid on
colorbar();

load('ind_field.mat')
figure(5), imagesc([0:360],[0:180],ind_field*1000)
ax = get(gca,'XTick');
set(gca,'XTick',[0:60:360])
set(gca,'GridLineStyle','-')
set(gca,'LineWidth',2)
grid on
colorbar();