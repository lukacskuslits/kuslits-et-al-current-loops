clear all, close all
raw_input = load('true_data_2010.mat');
raw_field = raw_input.true_data_2010.map;
raw_field = raw_field*1e-9;
raw_sv = raw_input.true_data_2010.sv;
raw_sv = raw_sv*1e-9;

% input = load('input_710.mat');
% raw_input = squeeze(input.vals(1,:,:));

orig_input = load('simulated_data.mat');
ref_input = orig_input.sim_data.maps_norm(:,:,122);

%index: 91 forrasok szama: 20 Merr: 0.09
%index: 10 forrasok szama: 35  Merr: 0.11
%index: 121 forrasok szama: 65 Merr: 0.17
%index: 710 forrasok szama: 85 Merr: 0.38
%find_label
% for sample = 1:1150
%     ref_input = orig_input.sim_data.maps_norm(:,:,sample);
%     xcorr = corr2(ref_input,raw_input);
%     if xcorr>0.7
%     disp(xcorr) 
%     disp(sample)
%     end
% end

output = load('output_2010.mat');
convert_back_to_dims('output_2010.mat')
outp_data = output.vals;

input = load('input_2010.mat');
inp_data = input.vals;
inp_data = squeeze(inp_data(1,:,:));
inp_sv = squeeze(input.vals(2,:,:));

outp_pos = squeeze(outp_data(1,1,:,:));
% 

% figure(2), imagesc(inp_data)
% figure(3), imagesc(inp_sv)
EstCoords=coordDerive(outp_pos,.2,.001,1);

% hold on,  scatter3(EstCoords(2,:),EstCoords(1,:),zed,100,'red','filled')
%EstCoords=coordDerive(outp_pos1.*abs(outp_I1),.2,.001,2);

eredmeny_control = load('eredmeny_control_oc_l_total_.mat');
% celfuggveny = load('celfuggveny_ls_total_.mat');
EstCoords = eredmeny_control.eredmeny_control(:,:,1,999);
% EstCoords(:,89) = [log(40);log(41)];
% EstCoords(:,14) = [log(8);log(80)];
% EstCoords(:,13) = [log(7);log(83)];
% EstCoords(:,16) = [log(9);log(63)];
% EstCoords(:,11) = [log(8);log(74)];
% EstCoords(:,10) = [log(5);log(66)];
% EstCoords(:,41) = [log(18);log(70)];
% EstCoords(:,21) = [log(12);log(65)];
% EstCoords(:,1) = [log(4);log(84)];
% EstCoords(:,6) = [log(4);log(62)];
[coef, pow] = transform_coefficients();
% %[res_est, sv_est] = EstFieldForGA(log(EstCoords), 45, 90, coef, pow);
% [res_static, sv_static] = EstFieldForGADims(log(EstCoords), 45, 90, 1, 0);
[res_raw, sv_raw] = EstFieldForGADims(EstCoords, 45, 90, 0, 0);
figure(2), imagesc(res_raw)

%% NRMS error 

ii = 1;
const_list = 1;
NRMS_LIST = zeros(length(const_list), 1);
MAE_LIST = zeros(length(const_list), 1);
xcorr_LIST = zeros(length(const_list), 1);

for constant = const_list 
    disp('ERTEKEK')
    disp(constant)
    [res_raw, sv_raw] = EstFieldForGADims(EstCoords, 45, 90, 0, constant);
    disp(max(max(res_raw)))
    %MAE = compute_NRMS_new(res,Val);
    [NRMS,MAE,xcorr] = computeNRMS(raw_sv,sv_raw);
    NRMS_LIST(ii) = NRMS; 
    MAE_LIST(ii) = MAE;
    xcorr_LIST(ii) = xcorr;
    ii = ii+1;
end

% NRMS = min(NRMS_LIST);
% constant = const_list(find(NRMS_LIST == NRMS));
%%    
% figure(4), imagesc(sv_raw)
% figure(5), imagesc(res_raw)
% figure(6), imagesc(sv_static)
% figure(7), imagesc(res_static)


dim_outp_pars = load('output_data_dims.mat');
dim_outp_pars = dim_outp_pars.output_data;
outp_depth = squeeze(dim_outp_pars(2,:,:));
outp_dtI = squeeze(dim_outp_pars(3,:,:));
outp_rad = squeeze(dim_outp_pars(4,:,:));
outp_I = squeeze(dim_outp_pars(5,:,:));


figure(5);
imagesc(outp_pos);
axesm ('mercator');
load coastlines;
plot(coastlon+180, coastlat+90, 'black');

figure(6);
imagesc(outp_depth);
axesm ('mercator');
load coastlines;
plot(coastlon+180, coastlat+90, 'black');

figure(7);
imagesc(outp_dtI);
clabel(inp_data_cont, k, 'LabelSpacing', 1000, 'FontSize', 6,'FontWeight','bold', 'Color', 'red')
axesm ('mercator');
load coastlines;
plot(coastlon+180, coastlat+90, 'black');

figure(8);
imagesc(outp_rad);
clabel(inp_data_cont, k, 'LabelSpacing', 1000, 'FontSize', 6,'FontWeight','bold', 'Color', 'red')
axesm ('mercator');
load coastlines;
plot(coastlon+180, coastlat+90, 'black');

figure(9);
imagesc(outp_I);
clabel(inp_data_cont, k, 'LabelSpacing', 1000, 'FontSize', 6,'FontWeight','bold', 'Color', 'red')
axesm ('mercator');
load coastlines;
plot(coastlon+180, coastlat+90, 'black');


%%
% figure(4), hold on 
% zed = 100*ones(84,1);
% scatter3(round(exp(EstCoords(2,:))),round(exp(EstCoords(1,:))),zed)
% figure(3), imagesc(res_raw)
% figure(4), imagesc(sv_raw)
% 
% save('est_field_cmb_final.mat','res_est')
% save('est_sv_cmb_final.mat','sv_est')
% save('est_field_cmb_raw.mat','res_raw')
% save('est_sv_cmb_raw.mat','sv_raw')
% 
% EstCoordsNew = exp(EstCoords);
% % zed = 10*ones(1,length(EstCoordsNew(1,:)));
% % figure(2), hold on,  scatter3(EstCoordsNew(2,:),EstCoordsNew(1,:),zed,100,'black','filled')
% 
% input = load('true_data.mat');%'input0.mat');
% inp_data = input.true_data.maps;
% inp_data1 = inp_data(:,:,1);%squeeze(inp_data(1,1,:,:));
% figure(5), imagesc(inp_data1)
% 
% sv_data = input.true_data.svs;
% sv_data1 = sv_data(:,:,1);%squeeze(inp_data(1,1,:,:));
% figure(6), imagesc(sv_data1)
% 
% save('true_field_cmb_final.mat','inp_data1')
% 
% inp_data_norm = input.true_data.maps_norm(:,:,1);
% sv_data_norm = input.true_data.svs_norm(:,:,1);
% sv_data_raw = input.true_data.svs(:,:,1);
% 
% 
% 
% %%
% %Quality metrics
% 
% disp(sum(sum(abs((res_raw-inp_data1))))/4050)
% disp(max(max(abs(inp_data1))))
% %max abs = 1.4e-3 T
% %mean abs 5e-4 T
% %MAE = 8.2e-4 T
% disp(sum(sum(abs((sv_raw-sv_data_raw))))/4050)
% disp(max(max(abs(sv_data_raw))))
% %max abs = 2.6e-13 T
% %mean abs 6.34e-14 T
% %MAE = 6.33e-14 T
% 
% 
% disp(corr2(inp_data_norm, res_est))
% disp(corr2(sv_data_norm, sv_est))
% %%
% % Surface field (only static!)
% % [res_est, res_est_raw] = EstFieldForGA(eredmeny_control.eredmeny_control(:,:,1,50), 45, 90);
% % figure(3), imagesc(res_est_raw)
