function convert_back_to_dimsPars(output_file)
dt = 25*3.15e7;
mu = 4*pi*1e-7;
attenuate = 1100;

rmin = 3.25;
rmax = 9;
% rmin = 0.25;%3.25; %3.25; %original upper and lower vals for R
% rmax = 1;%10; %10; %
Imin = 2.2e8;
Imax = 1e9;
%Imax = 5e7;%1e9;
depth_max = 3.48-0.15; %3.48e6;
depth_min = 2.7;
% depth_min = 2.7; %2.6; ICB
% depth_max = 3.3; %3.48; %3.48; CMB
dtB_max = 7.33e-13; %4e-11 with the lowest std
%Highest possible rate of change in current:
%dI/dt = 2*Rmax*(dB/dt)max/mu
dtI_max = 2*(((depth_max-depth_min)*1e6)^2+(3.25*1e5)^2)^(3/2)*(mu*(3.25*1e5)^2)^(-1)*dtB_max; %1.8182 [A/s]
dtI_max = dtI_max/attenuate; %Bind(dIt_max) ~ 5.59*Brmax

load(output_file)
if length(size(vals))==3
    depths = squeeze(vals(2,:,:));
    dtIs = squeeze(vals(3,:,:));
    rads = squeeze(vals(4,:,:));
    Is = squeeze(vals(5,:,:));
elseif length(size(vals))==4
    depths = squeeze(vals(1,2,:,:));
    dtIs = squeeze(vals(1,3,:,:));
    rads = squeeze(vals(1,4,:,:));
    Is = squeeze(vals(1,5,:,:));
end


depths = depths*(1e6*(depth_max-depth_min))+1e6*(depth_min-0.1);
depths(depths<1e6*(depth_min-0.1)) = 1e6*(depth_min-0.1);
depths(depths>1e6*depth_max) = 1e6*depth_max;
dtIs = dtIs*2*dtI_max - dtI_max;
rads = rads*1e5*rmax;
Is = Is*2*Imax-Imax;

output_data(2,:,:) = depths;
output_data(3,:,:) = dtIs;
output_data(4,:,:) = rads;
output_data(5,:,:) = Is;

save('output_data_dims.mat', 'output_data')
end