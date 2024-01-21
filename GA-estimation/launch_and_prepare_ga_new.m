clear all, close all
domx = 90;
domy = 180;
% To use for real data:
% convert_back_to_dims('output_2010.mat');
% input = load('input_2010.mat');
% output = load('output_2010.mat');
convert_back_to_dims('output_70.mat');
input = load('input_70.mat');
output = load('output_70.mat');
label = load('label_70.mat');
inp_data = input.vals;
outp_data = output.vals;
label_data = label.vals;
inp_data(1,:,:) = (inp_data(1,:,:) - min(min(inp_data(1,:,:))))/...
                  (max(max(inp_data(1,:,:))) - min(min(inp_data(1,:,:))));
inp_data(2,:,:) = (inp_data(2,:,:) - min(min(inp_data(2,:,:))))/...
                  (max(max(inp_data(2,:,:))) - min(min(inp_data(2,:,:))));
inp_data1 = squeeze(inp_data(1,:,:));
inp_data2 = squeeze(inp_data(2,:,:));
outp_pos1 = squeeze(outp_data(1,1,:,:));
outp_rad1 = squeeze(outp_data(1,2,:,:));
outp_I1 = squeeze(outp_data(1,3,:,:));
outp_depth1 = squeeze(outp_data(1,4,:,:));
outp_dtI1 = squeeze(outp_data(1,5,:,:));
%Compare with the true values (labels, only in case of synthetic data)
label_pos1 = squeeze(label_data(1,:,:));
label_rad1 = squeeze(label_data(2,:,:));
label_I1 = squeeze(label_data(3,:,:));
label_depth1 = squeeze(label_data(4,:,:));
label_dtI1 = squeeze(label_data(5,:,:));
f = inp_data;
save('ObjVal.mat','f')
% 
figure(1), 
contour(inp_data1)
% figure(2), 
% contour(inp_data2)
figure(2), 
subplot(1,2,1), imagesc(outp_pos1)
subplot(1,2,2), imagesc(label_pos1)

figure(3), 
subplot(1,2,1), imagesc(outp_rad1)
subplot(1,2,2), imagesc(label_rad1)

figure(4), 
subplot(1,2,1), imagesc(outp_I1)
subplot(1,2,2), imagesc(label_I1)

figure(5), 
subplot(1,2,1), imagesc(outp_depth1)
subplot(1,2,2), imagesc(label_depth1)

figure(6), 
subplot(1,2,1), imagesc(outp_dtI1)
subplot(1,2,2), imagesc(label_dtI1)


%EstCoords=coordDerive(outp_pos1.*abs(outp_I1),.2,.001,2);
deg_res = 2;
[coef, pow] = transform_coefficients_old(deg_res);
EstCoords=coordDerive(outp_pos1,.2,.001);

pars(1:10,:) = repmat(log(EstCoords),5,1);
pars(11:12,:) = 0.1;%randi(4605,2,length(EstCoords))/1000;
%pars(3,27) = 0.8;
[field_second, sv_second] = EstFieldForGANew(pars,domx,domy,coef,pow);
figure(8), contour(field_second)
% [field_first, sv_first] = EstFieldForGADims(log(EstCoords),domx,domy,0);
% figure(8), imagesc(field_first)
% display initial target for check
inp_data_norm = inp_data1;
inp_sv_norm = squeeze(inp_data(2,:,:));
target = sum(sum(abs(inp_data_norm-field_second)))/(domx*domy)+sum(sum(abs(inp_sv_norm-sv_second)))/(10*domx*domy);
disp(target)

EstCoords = repmat(EstCoords,5,1);
%prepare range arrays in ->/GA-eval
[l_array,h_array]=lh_array(EstCoords, 2, 90, 180);
zed = 100*ones(length(EstCoords(1,:)),1);
%figure(4), hold on, scatter3(EstCoords(2,:), EstCoords(1,:), zed, 100, 'black', 'filled')
figure(2), hold on, scatter3(l_array(2,:),l_array(1,:),zed,100,'red','filled')
figure(2), hold on, scatter3(h_array(2,:),h_array(1,:),zed,100,'blue','filled')

%prepare range
range=range_prepare(l_array,h_array);

%save range
range_log=log(range);
range_log(:,11) = repmat([0,6.9078],1,length(range_log(:,1))/2)';
range_log(:,12) = repmat([0,6.9078],1,length(range_log(:,1))/2)';
save('range_log.txt','range_log','-ascii')
% Don't forget to transform the range.txt with a header
%d a r R i I h H t T P p
%0 0 0 0 0 0 0 0 0 0 0 0
%And to separate the values in the file using only single spaces

close all, clear all
%% FINAL RECONSTRUCTION
subname=@EstFieldForGANew;
fparam=@params_loop;
[coef, pow] = transform_coefficients_old(2);
custom_genetic_algorithm(subname,fparam,'range_log.txt','ObjVal.mat',...
    12,3,8,32,90,180,coef,pow,'nloop',28,'mutrate',2/100,'d',0.13,...
    'const',1000,'migrate',1/8,'migtopo','ring','nii',501,'flmode','f',...
    'brcr',1e-3,'norma','rel','lconfig',16)