clear all, close all
output = load('output_70.mat');
outp_data = output.vals;
input = load('input_70.mat');
input_data = input.vals;
input_data(1,:,:) = (input_data(1,:,:) - min(min(input_data(1,:,:))))/...
                  (max(max(input_data(1,:,:))) - min(min(input_data(1,:,:))));
input_data(2,:,:) = (input_data(2,:,:) - min(min(input_data(2,:,:))))/...
                  (max(max(input_data(2,:,:))) - min(min(input_data(2,:,:))));
inp_data = squeeze(input_data(1,:,:));
inp_sv = squeeze(input_data(2,:,:));

% 
outp_pos = squeeze(outp_data(1,5,:,:));
figure(2), contour(inp_data, 50)
%figure(3), imagesc(inp_sv)
%EstCoords=coordDerive(outp_pos,.2,.001,1);

% hold on,  scatter3(EstCoords(2,:),EstCoords(1,:),zed,100,'red','filled')
%EstCoords=coordDerive(outp_pos1.*abs(outp_I1),.2,.001,2);

%Final GA result after the algorithm converges
% eredmeny_control = load('eredmeny_control_oc_l_total_.mat');
% EstCoords = eredmeny_control.eredmeny_control(:,:,500);
EstParCoords = load('final_ga_result.mat');
EstParCoords = EstParCoords.EstCoords;

[coef, pow] = transform_coefficients_old(2);
convert_back_to_dims('output_70.mat');
[res_est, sv_est] = EstFieldForGANew(EstParCoords, 90, 180, coef, pow);
%figure(4), imagesc(sv_est)
figure(5), contour(res_est, 50)
domx = 90;
domy = 180;
target = sum(sum(abs(inp_data-res_est)))/(domx*domy)+sum(sum(abs(inp_sv-sv_est)))/(10*domx*domy);
disp(target)
% compute static field - only in case of reconstructing the real
% geomagnetic map for total/static field rates
% [res_static, sv_static] = EstFieldForGADimsNew(EstCoords, 90, 180, 1);
% figure(6), imagesc(sv_static)
% figure(7), imagesc(res_static)
% figure(4), hold on 
% zed = 100*ones(84,1);
% scatter3(round(exp(EstCoords(2,:))),round(exp(EstCoords(1,:))),zed)

load('input_data_70_dims.mat');
load('output_data_dims.mat');
map_dims = true_data.maps;
for mult = 1.0:0.2:2.0
    disp('MULT: ')
    disp(mult)
    tmp = output_data(3,:,:);
    output_data(3,:,:) = mult*tmp;
    save('output_data_dims.mat', 'output_data');
    [res_raw, sv_raw] = EstFieldForGADimsNew(EstParCoords, 90, 180, 0);
    disp(sum(sum(abs((map_dims-res_raw))))/16400)
end
convert_back_to_dims('output_70.mat');
load('output_data_dims.mat');
output_data(3,:,:) = 2.8*output_data(3,:,:);
save('output_data_dims.mat', 'output_data');

[res_raw, sv_raw] = EstFieldForGADimsNew(EstParCoords, 90, 180, 0);
% res_raw = interpn(res_raw, 1);
% sv_raw = interpn(sv_raw, 1);
% figure(8), contour(res_raw, 50)
% figure(9), contour(sv_raw, 50)


save('est_field_cmb_synth.mat','res_est')
save('est_sv_cmb_synth.mat','sv_est')
save('est_field_cmb_synth_raw.mat','res_raw')
save('est_sv_cmb_synth_raw.mat','sv_raw')

%EstCoordsNew = exp(EstParCoords);
% zed = 10*ones(1,length(EstCoordsNew(1,:)));
% figure(2), hold on,  scatter3(EstCoordsNew(2,:),EstCoordsNew(1,:),zed,100,'black','filled')

input = load('input_data_70_dims.mat');%'input0.mat');
inp_data = input.true_data.maps;
inp_data1 = inp_data(:,:,1);%squeeze(inp_data(1,1,:,:));
figure(5), imagesc(inp_data1)

sv_data = input.true_data.svs;
sv_data_raw = sv_data(:,:,1);%squeeze(inp_data(1,1,:,:));
figure(6), imagesc(sv_data_raw)

inp_data_norm = inp_data;
sv_data_norm = inp_sv;

%%
%%Quality metrics
% parameter error (in case of synthetic data)
inp_data_label = load('label_70.mat');
convert_back_to_dims('label_70.mat');
inp_pos = squeeze(inp_data_label.vals(1, :, :));
inp_Is = squeeze(inp_data_label.vals(3, :, :));
figure(13), imagesc(inp_Is);
EstCoords=coordDerive(inp_pos,.2,.001);
EstCoords(:,34) = [1;9];
EstCoords(:,35) = [1;82];
TruePars = get_true_pars(EstCoords);

EstPars = get_est_pars(EstParCoords, 90, 180);
[merr, prefc] = calculate_parameter_error_loop(EstPars, TruePars, 8, 28);
Merr = mean(mean(merr));

% misfit values
disp(sum(sum(abs((map_dims-res_raw))))/16400)
disp(max(max(abs(res_raw))))
disp(max(max(abs(map_dims))))

%max abs =  7.3e-4 T
%MAE = 1.106e-4 T; 
%NRMS = 0.054
disp(sum(sum(abs((sv_raw-sv_data_raw))))/16400)
disp(max(max(abs(sv_data_raw))))
%max abs = 2.6e-13 T
%mean abs 6.34e-14 T
%MAE = 6.33e-14 T


disp(corr2(inp_data, res_est))
disp(corr2(inp_sv, sv_est))
%%
% Surface field (only static!)
% [res_est, res_est_raw] = EstFieldForGA(eredmeny_control.eredmeny_control(:,:,1,50), 45, 90);
% figure(3), imagesc(res_est_raw)
