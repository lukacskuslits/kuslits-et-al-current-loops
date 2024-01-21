function [EstField, sv_i] = EstFieldForGADimsNew(pars, domx, domy, static)
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

EstParsOrig = EstPars;

%disp(EstPars)
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
