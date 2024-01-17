function [sim_data, sample_numbers] = test_train_generator_new(attenuate, deg_res)
%% Generates a training set with a specified number of current loops and batch size
%Input params:
% attenuate: 
% deg_res:
% set: number of loops in the training set batches
% batch_size: size of a training set batch

%orig outlp arguments: [sim_data, true_data, sample_batch_sizes]
%clear all, close all;
%    mu = 4*pi*1e-7;
%     rmax = 9e5;
%     rmin = 3.25e5;
%     Imax = 1e9;
%     depth_max = 3.48e6-1.5e5; %3.48e6;
%     depth_min = 2.7e6;
%     dBt_max = 7.3e-13; %4e-11 with the lowest std
%     %Highest possible rate of change in current:
%     %dI/dt = 2*Rmax*(dB/dt)max/mu
%     dIt_max = 2*((depth_max-depth_min)^2+rmin^2)^(3/2)*(mu*rmin^2)^(-1)*dBt_max;
%     dIt_max = dIt_max/attenuate;
    lat_dim = 180/deg_res;
    long_dim = 360/deg_res;
    % maps of the field and secular variation values
    simulated_field_maps_raw = zeros(lat_dim, long_dim, 1);
    simulated_field_maps_norm = zeros(lat_dim, long_dim, 1);
    simulated_field_maps_log = zeros(lat_dim, long_dim, 1);
    simulated_sv_maps_raw = zeros(lat_dim, long_dim, 1);
    simulated_sv_maps_norm = zeros(lat_dim, long_dim, 1);
    
    % maps of source positions and source parameters
    simulated_pos_data_total = zeros(lat_dim, long_dim, 1);
    simulated_rad_maps_total = zeros(lat_dim, long_dim, 1);
    simulated_I_maps_total = zeros(lat_dim, long_dim, 1);
    simulated_depth_maps_total = zeros(lat_dim, long_dim, 1);
    simulated_dtI_maps_total = zeros(lat_dim, long_dim, 1);
%     simulated_vlat_maps_total = zeros(lat_dim, long_dim, 1);
%     simulated_vlong_maps_total = zeros(lat_dim, long_dim, 1);
    
    %set = [25, 60, 50, 40, 50, 100, 400, 1000];
    %set = 25:100:1025;
    %     y = makedist('Uniform', 'lower', 25, 'upper', 1025);
    %     y = pdf(y,x);
    %     set = 25:10:50;
    %     x = 1:50;
    %     y = normpdf(x, 35, 10);
%     set = 5:5:550;
%     x = 1:550;
%     %training set generated
%     %lat_dim, 2
%     %65, 2
%     %110, 2
%     %120, 2
%     y = normpdf(x, 65, 2); 
%     batch_size = round(2e2*y(set)); %default to: 1e3
%     disp(batch_size)
%     set = set(batch_size~=0);
%     batch_size = batch_size(batch_size~=0);
%     disp('Set:')
%     disp(set)
%     disp('batch_size:')
%     disp(batch_size)
    set = 1:5:6;%25:10:175; %;%55:10:155;
    batch_size = 2;
    sample_numbers = zeros(1,length(set));
    for ii = 1:length(set)
        disp('Nr. of sources in set:')
        disp(set(ii))
        try
            [pos_data_all, rad_data_all, I_data_all, depth_data_all, dtI_data_all, br_all, sv_all, rnorm_all, sv_norm_all] = gen_loops_par_rad_total_new(batch_size, set(ii), attenuate, deg_res);
        catch ME
           disp(ME)
           break
        end
        sample_number = 0;
        for sample = 1:length(br_all(1,1,:))
%              disp('Index of inspected sample:')
%              disp(sample)
            % select only samples which are inside the fences Q1 - 3*IQ
            % and Q3 + 3*IQ determined using true field model data 
%             disp('Bmax')
%             disp(max(max(br_all(:,:,sample))))
%             disp('Bmin')
%             disp(min(min(br_all(:,:,sample))))
            if (max(max(br_all(:,:,sample)))<=0.0032) && (min(min(br_all(:,:,sample)))>=-0.0031)
                %disp('Sample accepted!')
                br = br_all(:,:,sample);
                pos_data = pos_data_all(:,:,sample);
                rad_data = rad_data_all(:,:,sample);
                depth_data = depth_data_all(:,:,sample);
                I_data = I_data_all(:,:,sample);
                dtI_data = dtI_data_all(:,:,sample);
                sv = sv_all(:,:,sample);
                sv_norm = sv_norm_all(:,:,sample);
                rnorm = rnorm_all(:,:,sample);
                simulated_field_maps_raw = cat(3, simulated_field_maps_raw, br);
                simulated_field_maps_norm = cat(3, simulated_field_maps_norm, rnorm);
                simulated_field_maps_log = cat(3, simulated_field_maps_log, log(abs(rnorm)));
                simulated_sv_maps_raw = cat(3, simulated_sv_maps_raw, sv);
                simulated_sv_maps_norm = cat(3, simulated_sv_maps_norm, sv_norm);
                simulated_pos_data_total = cat(3, simulated_pos_data_total, pos_data);
                simulated_rad_maps_total = cat(3, simulated_rad_maps_total, rad_data);
                simulated_I_maps_total = cat(3, simulated_I_maps_total, I_data);
                simulated_depth_maps_total = cat(3, simulated_depth_maps_total, depth_data);
                simulated_dtI_maps_total = cat(3, simulated_dtI_maps_total, dtI_data);
                sample_number = sample_number + 1;
            end
        end
        sample_numbers(ii) = sample_number;
        if sample_number == 0
            disp('No sample could be accepted')
            break
        end  
%         simulated_vlat_maps_total = cat(3, simulated_vlat_maps_total, vlat_data);
%         simulated_vlong_maps_total = cat(3, simulated_vlong_maps_total, vlong_data);
    end
    
    simulated_pos_data_total = simulated_pos_data_total(:,:,2:end);
    simulated_field_maps_raw = simulated_field_maps_raw(:,:,2:end);
    simulated_field_maps_norm = simulated_field_maps_norm(:,:,2:end);
    simulated_field_maps_log = simulated_field_maps_log(:,:,2:end);
    simulated_sv_maps_raw = simulated_sv_maps_raw(:,:,2:end);
    simulated_sv_maps_norm = simulated_sv_maps_norm(:,:,2:end);
    simulated_rad_maps_total = simulated_rad_maps_total(:,:,2:end);
    simulated_I_maps_total = simulated_I_maps_total(:,:,2:end);
    simulated_depth_maps_total = simulated_depth_maps_total(:,:,2:end);
    simulated_dtI_maps_total = simulated_dtI_maps_total(:,:,2:end);
%     simulated_vlat_maps_total = simulated_vlat_maps_total(:,:,2:end);
%     simulated_vlong_maps_total = simulated_vlong_maps_total(:,:,2:end);

    
%% Estimate field values and SV for forecasting
    dt = 5*3.15e7;
    true_field_maps_raw = zeros(lat_dim, long_dim, 1);
    true_field_maps_norm = zeros(lat_dim, long_dim, 1);
    true_field_maps_log = zeros(lat_dim, long_dim, 1);
    
    test = load('test_core_new.txt','-ascii');
    ii = 1;
    for time = 1600:5:1975
        disp(time)
        test_i = test(test(:,1)==time,:);
        if time==1815
            test_i=test_i(1:lat_dim*long_dim,:);
        end
        BR = test_i(:,4);
        length(BR)
        BR_norm = (BR-min(BR))./(max(BR)-min(BR));
        BR = reshape(-BR,long_dim, lat_dim)';
        BR = shiftarray(BR, lat_dim, 'col');
        BR = mirror(BR,'row','surfa');
        BR_norm = reshape(-BR_norm,long_dim, lat_dim)';
        BR_norm = shiftarray(BR_norm, lat_dim, 'col');
        BR_norm = mirror(BR_norm,'row','surfa');
        true_field_maps_raw = cat(3, true_field_maps_raw, BR/1e9);
        true_field_maps_norm = cat(3, true_field_maps_norm, BR_norm);
        true_field_maps_log = cat(3, true_field_maps_log, log(abs(BR_norm)));
        ii = ii+1;
    end    
    
true_field_maps_raw = true_field_maps_raw(:,:,2:end);
true_field_maps_norm = true_field_maps_norm(:,:,2:end);
true_field_maps_log = true_field_maps_log(:,:,2:end);
diffmap = diff(true_field_maps_raw,1,3);
true_sv_maps_raw = diffmap/dt;
min_svs = min(min(true_sv_maps_raw));
max_svs = max(max(true_sv_maps_raw));
true_sv_maps_norm = bsxfun(@rdivide,bsxfun(@minus,true_sv_maps_raw,min_svs),(max_svs-min_svs));
true_field_maps_raw = true_field_maps_raw(:,:,2:end);
true_field_maps_norm = true_field_maps_norm(:,:,2:end);
true_field_maps_log = true_field_maps_log(:,:,2:end);
    
%randomly shuffle training set
ssim = size(simulated_pos_data_total);
idx = randperm(ssim(3));
sim_pos = simulated_pos_data_total;
sim_pos(:,:,idx) = simulated_pos_data_total(:,:,:);
sim_map = simulated_field_maps_raw;
sim_map(:,:,idx) = simulated_field_maps_raw(:,:,:);
sim_normmap = simulated_field_maps_norm;
sim_normmap(:,:,idx) = simulated_field_maps_norm(:,:,:);
sim_logmap = simulated_field_maps_log;
sim_logmap(:,:,idx) = simulated_field_maps_log(:,:,:);
sim_rad = simulated_rad_maps_total;
sim_rad(:,:,idx) = simulated_rad_maps_total(:,:,:);
sim_I = simulated_I_maps_total;
sim_I(:,:,idx) = simulated_I_maps_total(:,:,:);
sim_depth = simulated_depth_maps_total;
sim_depth(:,:,idx) = simulated_depth_maps_total(:,:,:);
sim_dtI = simulated_dtI_maps_total;
sim_dtI(:,:,idx) = simulated_dtI_maps_total(:,:,:);
% sim_vlat = simulated_vlat_maps_total;
% sim_vlat(:,:,idx) = simulated_vlat_maps_total(:,:,:);
% sim_vlong = simulated_vlong_maps_total;
% sim_vlong(:,:,idx) = simulated_vlong_maps_total(:,:,:);
sim_sv = simulated_sv_maps_raw;
sim_sv(:,:,idx) = simulated_sv_maps_raw(:,:,:);
sim_svnorm = simulated_sv_maps_norm;
sim_svnorm(:,:,idx) = simulated_sv_maps_norm(:,:,:);
labels_sim = zeros(1,ssim(3));


%prepare labelled data sets as structs
strue = size(true_field_maps_log);
idx = randperm(strue(3));
true_map = true_field_maps_raw;
true_map(:,:,idx) = true_field_maps_raw(:,:,:);
true_normmap = true_field_maps_norm;
true_normmap(:,:,idx) = true_field_maps_norm(:,:,:);
true_logmap = true_field_maps_log;
true_logmap(:,:,idx) = true_field_maps_log(:,:,:);
true_sv = true_sv_maps_raw;
true_sv(:,:,idx) = true_sv_maps_raw(:,:,:);
true_sv_norm = true_sv_maps_norm;
true_sv_norm(:,:,idx) = true_sv_maps_norm(:,:,:);
labels_true = ones(1,strue(3));
true_data.labels = labels_true;
true_data.maps = true_map;
true_data.maps_norm = true_normmap;
true_data.logmaps = true_logmap;
true_data.svs = true_sv;
true_data.svs_norm = true_sv_norm;

sim_data.labels = labels_sim;
sim_data.maps = sim_map;
sim_data.maps_norm = sim_normmap;
sim_data.logmaps = sim_logmap;
sim_data.svs = sim_sv;
sim_data.svs_norm = sim_svnorm;
sim_data.pos = sim_pos;
sim_data.rads = sim_rad;
sim_data.Is = sim_I;
sim_data.depths = sim_depth;
sim_data.dtIs = sim_dtI;
% sim_data.vlats = sim_vlat;
% sim_data.vlongs = sim_vlong;

% save('simulated_data_new.mat', 'sim_data');
% save('true_data_new.mat','true_data');
% save('accepted_samples_new.mat','sample_batch_sizes')
end
