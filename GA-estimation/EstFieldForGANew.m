function [EstFieldNorm, sv_norm] = EstFieldForGANew(pars, domx, domy, coef, pow)
pars(11:12,:) = -pars(11:12,:);
pars = exp(pars);

EstCoords = round(pars(1:10,:));
chance_arrI = pars(11,:);
chance_arrDtI = pars(12,:);

[pcal_lat, pcal_lon] = calcCoordBounds(EstCoords, domx, domy);
EstPosCoords(1,:) = pcal_lat(1,:);
EstPosCoords(2,:) = pcal_lon(1,:);
EstRadCoords(1,:) = pcal_lat(2,:);
EstRadCoords(2,:) = pcal_lon(2,:);
EstICoords(1,:) = pcal_lat(3,:);
EstICoords(2,:) = pcal_lon(3,:);
EstDepthCoords(1,:) = pcal_lat(4,:);
EstDepthCoords(2,:) = pcal_lon(4,:);
EstDtICoords(1,:) = pcal_lat(5,:);
EstDtICoords(2,:) = pcal_lon(5,:);

outp_data = load('output_data_dims.mat');
outp_data = outp_data.output_data;
outp_rad1 = squeeze(outp_data(2,:,:));
outp_I1 = squeeze(outp_data(3,:,:));
outp_depth1 = squeeze(outp_data(4,:,:));
outp_dtI1 = squeeze(outp_data(5,:,:));

EstRadCoords = sub2ind(size(outp_depth1),EstRadCoords(1,:),EstRadCoords(2,:));
EstICoords = sub2ind(size(outp_depth1),EstICoords(1,:),EstICoords(2,:));
EstDepthCoords = sub2ind(size(outp_depth1),EstDepthCoords(1,:),EstDepthCoords(2,:));
EstDtICoords = sub2ind(size(outp_depth1),EstDtICoords(1,:),EstDtICoords(2,:));

EstRads = outp_rad1(EstRadCoords);
EstIs = outp_I1(EstICoords);
EstDepths = outp_depth1(EstDepthCoords);
EstDtIs = outp_dtI1(EstDtICoords);


%current and current change rate flip by chance
EstIs(chance_arrI>0.5) = -EstIs(chance_arrI>0.5);
EstDtIs(chance_arrDtI>0.5) = -EstDtIs(chance_arrDtI>0.5);

EstPars = [EstPosCoords;EstRads;EstIs;EstPosCoords;EstDepths;EstDtIs];
EstPars(5,:)=2*EstPars(1,:)*pi()/180;
EstPars(6,:)=2*EstPars(2,:)*pi()/180;
EstPars(1,:)=2*EstPars(1,:)*pi()/180;
EstPars(2,:)=2*EstPars(2,:)*pi()/180;

disp(EstPars)
EstField = total_field_approx_new_ga(EstPars, domx, domy, coef, pow);

EstField = reshape(EstField, domy, domx)';
EstFieldNorm = (EstField - min(min(EstField)))/(max(max(EstField))-min(min(EstField)));

EstParsOrig = EstPars;
EstParsOrig(8,:)=0;
dt = 5*3.15e7;
EstFieldOrig = total_field_approx_new_ga(EstParsOrig, domx, domy, coef, pow);
EstFieldOrig = reshape(EstFieldOrig, domy, domx)';
pars_new = EstParsOrig;
pars_new(4,:) = pars_new(4,:) + EstPars(8,:)*dt;
res_new = total_field_approx_new_ga(pars_new, domx, domy, coef, pow);
res_new = reshape(res_new, domy, domx)';
sv_i = (res_new - EstFieldOrig)/dt;
sv_norm = (sv_i-min(min(sv_i)))./(max(max(sv_i))-min(min(sv_i)));
end
% 
% MSV = std(std(sv_i));
% for ii=1:45
%     for jj = 1:90
%         if abs(sv_i(ii,jj))>8*MSV
%            sv_i(ii,jj) = MSV;
%         end
%     end
% end
% EstPars(4,11) = -2e9;
% EstPars(3,11) = 8e5;
% EstPars(7,11) = 3.8e6;
% EstPars(8,11) = -1e-14;
% EstPars(4,21) = -2e9;
% EstPars(3,21) = 8e5;
% EstPars(7,21) = 3.8e6;
% EstPars(8,21) = -1e-14;
% EstPars(4,14) = -1e9;
% EstPars(3,14) = 8e5;
% EstPars(7,14) = 3.7e6;
% EstPars(8,14) = 1e-14;
% EstPars(4,16) = -1e9;
% EstPars(3,16) = 8e5;
% EstPars(7,16) = 3.7e6;
% EstPars(8,16) = 1e-14;
% EstPars(4,13) = -1e9;
% EstPars(3,13) = 8e5;
% EstPars(7,13) = 3.7e6;
% EstPars(8,13) = 5e-14;
% EstPars(4,1) = -1e9;
% EstPars(3,1) = 8e5;
% EstPars(7,1) = 3.7e6;
% EstPars(8,1) = -1e-14;
% EstPars(4,10) = -3e9;
% EstPars(3,10) = 8e5;
% EstPars(7,10) = 3.8e6;
% EstPars(8,10) = 1e-14;
% 
% 
% EstPars(4,72) = 2.5e9;
% EstPars(3,72) = 1e6;
% EstPars(7,72) = 3.8e6;
% EstPars(8,72) = -1e-14;
% EstPars(4,78) = 2e9;
% EstPars(3,78) = 1e6;
% EstPars(7,78) = 3.8e6;
% EstPars(8,78) = -1e-14;
% EstPars(4,79) = 1e9;
% EstPars(3,79) = 1.5e6;
% EstPars(7,79) = 3.8e6;
% EstPars(8,79) = -1e-14;
% EstPars(4,90) = 1e9;
% EstPars(3,90) = 8e5;
% EstPars(7,90) = 3.8e6;
% EstPars(8,90) = -1e-14;
% EstPars(4,88) = 4e8;
% EstPars(3,88) = 1e6;
% EstPars(7,88) = 3.8e6;
% EstPars(8,88) = -1e-14;
% EstPars(4,84) = 4e8;
% EstPars(3,84) = 1e6;
% EstPars(7,84) = 3.8e6;
% EstPars(8,84) = -1e-14;
% EstPars(4,87) = 1e9;
% EstPars(3,87) = 1e6;
% EstPars(7,87) = 3.8e6;
% EstPars(8,87) = -1e-14;
% 
