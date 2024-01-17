function convert_back_to_dims(output_file)
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



outp_data = load(output_file);
if length(size(outp_data))==4
    rads = squeeze(outp_data.vals(1,2,:,:));
    Is = squeeze(outp_data.vals(1,3,:,:));
    depths = squeeze(outp_data.vals(1,4,:,:));
    dtIs = squeeze(outp_data.vals(1,5,:,:));
else
    rads = squeeze(outp_data.vals(2,:,:));
    Is = squeeze(outp_data.vals(3,:,:));
    depths = squeeze(outp_data.vals(4,:,:));
    dtIs = squeeze(outp_data.vals(5,:,:));
end

rads(rads<0)=0;
depths(depths<0)=0;

    
depths = depths*(1e6*(depth_max-depth_min))+1e6*(depth_min);
depths(depths<1e6*(depth_min)) = 1e6*(depth_min);
depths(depths>1e6*depth_max) = 1e6*depth_max;
dtIs = dtIs*2*dtI_max - dtI_max;
rads = rads*1e5*rmax;
Is = Is*2*Imax-Imax;
disp('Rads_d')
disp(min(min(rads)))

output_data(2,:,:) = rads;
output_data(3,:,:) = Is;
output_data(4,:,:) = depths;
output_data(5,:,:) = dtIs;

save('output_data_dims.mat', 'output_data')
end