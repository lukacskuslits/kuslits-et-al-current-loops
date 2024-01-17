function [EstFieldNorm, EstPars] = EstFieldForGAPars(EstCoords, domx, domy, coef, pow)
EstCoords = round(exp(EstCoords));
disp(EstCoords)
if EstCoords<0.5
   EstCoords=0.5;
end
pcal_lat = EstCoords(1,:);
pcal_lon = EstCoords(2,:);
pcal_lat(pcal_lat>domx) = domx;
pcal_lon(pcal_lon>domy) = domy;

EstCoords(1,:) = pcal_lat;
EstCoords(2,:) = pcal_lon;
% disp('MIERT ROSSZ?')
% disp(EstCoords)

outp_data = load('output_data_dims.mat');
outp_data = outp_data.output_data;
disp(size(outp_data))
outp_depth1 = squeeze(outp_data(2,:,:));
outp_dtI1 = squeeze(outp_data(3,:,:));
outp_rad1 = squeeze(outp_data(4,:,:));
outp_I1 = squeeze(outp_data(5,:,:));
RefCoords = sub2ind(size(outp_depth1),EstCoords(1,:),EstCoords(2,:));

EstDepths = outp_depth1(RefCoords);
EstDtIs = outp_dtI1(RefCoords);
EstRads = outp_rad1(RefCoords);
EstIs = outp_I1(RefCoords);

EstPars = [EstCoords;EstRads;EstIs;EstCoords;EstDepths;EstDtIs];

EstPars(5,:)=4*EstPars(1,:)*pi()/180;
EstPars(6,:)=4*EstPars(2,:)*pi()/180;
EstPars(1,:)=4*EstPars(1,:)*pi()/180;
EstPars(2,:)=4*EstPars(2,:)*pi()/180;

EstField = total_field_approx(EstPars, domx, domy, coef, pow);

EstField = reshape(EstField, domy, domx)';
EstFieldNorm = (EstField - min(min(EstField)))/(max(max(EstField))-min(min(EstField)));

EstParsOrig = EstPars;
EstParsOrig(8,:)=0;
dt = 5*3.15e7;
EstFieldOrig = total_field_approx(EstParsOrig, domx, domy, coef, pow);
EstFieldOrig = reshape(EstFieldOrig, domy, domx)';
pars_new = EstParsOrig;
pars_new(4,:) = pars_new(4,:) + EstPars(8,:)*dt;
res_new = total_field_approx(pars_new, domx, domy, coef, pow);
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
