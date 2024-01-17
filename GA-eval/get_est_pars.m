function [EstPars] = get_est_pars(pars, domx, domy)
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

