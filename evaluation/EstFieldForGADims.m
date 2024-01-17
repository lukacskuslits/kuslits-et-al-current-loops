function [EstField, sv_i] = EstFieldForGADims(EstCoords, domx, domy, static)
EstCoords = round(exp(EstCoords));
if EstCoords<0.5
   EstCoords=0.5;
end
pcal_lat = EstCoords(1,:);
pcal_lon = EstCoords(2,:);
pcal_lat(pcal_lat>domy) = domy;
pcal_lon(pcal_lon>domx) = domx;

EstCoords(1,:) = pcal_lat;
EstCoords(2,:) = pcal_lon;

outp_data = load('output_data_dims.mat');
outp_data = outp_data.output_data;
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

EstPars(5,:)=2*EstPars(1,:)*pi()/180;
EstPars(6,:)=2*EstPars(2,:)*pi()/180;
EstPars(1,:)=2*EstPars(1,:)*pi()/180;
EstPars(2,:)=2*EstPars(2,:)*pi()/180;

EstParsOrig = EstPars;

disp(EstPars)
if static==1
    EstPars(8,:)=0;
    EstParsOrig(8,:)=0;
end
EstField = total_field_approx_new(EstPars, domx, domy);
EstField = reshape(EstField, domy, domx)';
%EstFieldNorm = (EstField - min(min(EstField)))/(max(max(EstField))-min(min(EstField)));

dt = 5*3.15e7;
EstFieldOrig = total_field_approx_new(EstParsOrig, domx, domy);
EstFieldOrig = reshape(EstFieldOrig, domy, domx)';
pars_new = EstParsOrig;
pars_new(4,:) = pars_new(4,:) + EstParsOrig(8,:)*dt;
res_new = total_field_approx_new(pars_new, domx, domy);
res_new = reshape(res_new, domy, domx)';
sv_i = (res_new - EstFieldOrig)/dt;
%MSV = std(std(sv_i));
% for ii=1:45
%     for jj = 1:90
%         if abs(sv_i(ii,jj))>8*MSV
%            sv_i(ii,jj) = MSV;
%         end
%     end
% end
%sv_norm = (sv_i-min(min(sv_i)))./(max(max(sv_i))-min(min(sv_i)));
end
% 
