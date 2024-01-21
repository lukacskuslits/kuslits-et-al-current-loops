function [pos_data, rad_data, I_data, depth_data, dtI_data, res, sv, rnorm, sv_norm] = gen_loops_par_rad_total_new(set_size,nloops,attenuate, deg_res)
%% Generates a training set containing a varying set_size of current loops placed inside the outer core
%% with differing loop parameters assigned randomly within the parameter ranges specified using REF!.m
%Input params: 
%-------------
% set_size: set_size of loop models in the given batch of the training set
% nloops: set_size of loops in each model in the batch
% attenuate: factor (\gamma) for reducing the maximum potential rate of change in the loop currents (dI/dt)
% deg_res: resolution of output field maps in spherical coordinates (degrees)
%Input files:
%------------
% coef_new.mat: matrix containing polynomial coeffitients used for the
    % approximations of diffusion screening
% powers_new.mat: matrix containing polynomial powers used for the
    % approximations of diffusion screening
%Output values:
%--------------
% pos_data: distribution of loop positions
% rad_data: distribution of loop radii
% I_data: distribution of loop current intensities
% depth_data: distribution of loop depths
% dtI_data: distribution of dI/dt values
% res: radial field data
% sv: radial SV data
% rnorm: normalized radial field data
% sv_norm: normalized radial SV data


dt = 25*3.15e7;
mu = 4*pi*1e-7;
lat_dim = 180/deg_res;
long_dim = 360/deg_res;

transform_coefficients_new(deg_res);

pos_data=zeros(lat_dim, long_dim, set_size);
rad_data=zeros(lat_dim, long_dim, set_size);
I_data=zeros(lat_dim, long_dim, set_size);
depth_data=zeros(lat_dim, long_dim, set_size);
dtI_data=zeros(lat_dim, long_dim, set_size);
res=zeros(lat_dim, long_dim,set_size);
sv=zeros(lat_dim, long_dim, set_size);
rnorm=zeros(lat_dim, long_dim, set_size);
sv_norm=zeros(lat_dim, long_dim, set_size);


rmin = 3.25;
rmax = 9;
Imin = 2.2e8;
Imax = 1e9;
depth_max = 3.48-0.15; %3.48e6;
depth_min = 2.7;
dtB_max = 7.33e-13; %4e-11 with the lowest std
%Highest possible rate of change in current:
dtI_max = 2*(((depth_max-depth_min)*1e6)^2+(3.25*1e5)^2)^(3/2)*(mu*(3.25*1e5)^2)^(-1)*dtB_max; %1.8182 [A/s]
dtI_max = dtI_max/attenuate; %Bind(dIt_max) ~ 5.59*Brmax
dtI_min = 0;
v_max = 0;

%
parfor ii=1:set_size
    limits = [1;lat_dim-1;1;long_dim-1;rmin;rmax;Imin;Imax;-1;1;0;0;0;0;depth_min;depth_max;dtI_min;dtI_max;-v_max;v_max];
    try
        pp = loop_gen_constr_rad_total(limits,nloops,3,3,3,3,[deg_res,deg_res]);     
    catch ME
        disp(ME)
        res(:,:,ii) = ones(lat_dim, long_dim);
        continue
    end
    try 
        % generate source detection probability mask with simulated gauss distribution for ML
        cords=pp(1:2,:)/deg_res;
        [pmarks,radmarks,Imarks_neg,Imarks_pos,depthmarks,dtImarks_neg,dtImarks_pos,vlats_neg,vlats_pos,vlongs_neg,vlongs_pos]=pos_marks([cords;pp(3:end,:)],lat_dim,long_dim,1,1);
        v = [5, 5];
        h = [5, 5];
        pos_data_i = gaussing(pmarks, 3, v, h);
        rad_data_i = gaussing(radmarks, 3, v, h);
        I_data_i_neg = gaussing(Imarks_neg, 3, v, h);
        I_data_i_pos = gaussing(Imarks_pos, 3, v, h);
        I_data_i = I_data_i_neg + I_data_i_pos;
        depth_data_i = gaussing(depthmarks, 3, v, h);
        dtI_data_i_neg = gaussing(dtImarks_neg, 3, v, h);
        dtI_data_i_pos = gaussing(dtImarks_pos, 3, v, h);
        dtI_data_i = dtI_data_i_neg + dtI_data_i_pos;
%         vlat_data_i_neg = gaussing(vlats_neg, 3, v, h);
%         vlat_data_i_pos = gaussing(vlats_pos, 3, v, h);
%        vlat_data_i = vlat_data_i_neg + vlat_data_i_pos;
%         vlong_data_i_neg = gaussing(vlongs_neg, 3, v, h);
%         vlong_data_i_pos = gaussing(vlongs_pos, 3, v, h);
%        vlong_data_i = vlong_data_i_neg + vlong_data_i_pos;
        pos_data(:,:,ii) = pos_data_i;
        rad_data(:,:,ii) = rad_data_i/(1e5*rmax);
        I_data(:,:,ii) = (I_data_i+Imax)/(2*Imax);
        depth_data(:,:,ii) = (depth_data_i-1e6*depth_min)/(1e6*(depth_max-depth_min));
        dtI_data(:,:,ii) = (dtI_data_i + dtI_max)/(2*dtI_max);
%         vlat_data(:,:,ii) = vlat_data_i/v_max;
%         vlong_data(:,:,ii) = vlong_data_i/v_max;

        %normalized radial field at CMB depth
        %generated result spherical coordinates
        pp(5,:)=pp(1,:)*pi()/180;
        pp(6,:)=pp(2,:)*pi()/180;
        pp(1,:)=pp(1,:)*pi()/180;
        pp(2,:)=pp(2,:)*pi()/180;
        res_i = total_field_approx_new(pp, lat_dim, long_dim);
        rnorm_i = (res_i-min(res_i))./(max(res_i)-min(res_i));
    %   rnorm(:,ii)=(rnorm1-min(rnorm1))./(max(rnorm1)-min(rnorm1));
        rnorm_i = reshape(rnorm_i, long_dim, lat_dim);
        rnorm_i = rnorm_i';
        rnorm(:,:,ii) = rnorm_i;
        res(:,:,ii) = reshape(res_i, long_dim, lat_dim)';

        %secular variation of the radial field during time dt
        pars_new = pp;
        pars_new(4,:) = pars_new(4,:) + pp(8,:)*dt;
        pars_new(1,:) = pars_new(1,:) + pp(9,:)*dt;
        pars_new(2,:) = pars_new(2,:) + pp(10,:)*dt;
        cord_max = [lat_dim,long_dim];
        for cord_part=1:2
            pars_part = pars_new(cord_part,:);
            pars_part(pars_part<0) = cord_max(cord_part)-pars_part(pars_part<0);
            pars_new(cord_part,:) = pars_part;
        end
        res_new = total_field_approx_new(pars_new, lat_dim, long_dim);
        sv_i = (res_new - res_i)/dt;
        sv(:,:,ii) = reshape(sv_i, long_dim, lat_dim)';
        svnorm_i = (sv_i-min(sv_i))./(max(sv_i)-min(sv_i));
        sv_norm(:,:,ii) = reshape(svnorm_i, long_dim, lat_dim)';
    catch ME
        disp(getReport(ME,'extended'))
        rethrow(ME)
        poolobj = gcp('nocreate');
        delete(poolobj)
    end
end

end
