function [TruePars] = get_true_pars(EstCoords)

outp_data = load('output_data_dims.mat');
outp_data = outp_data.output_data;
outp_rad1 = squeeze(outp_data(2,:,:));
outp_I1 = squeeze(outp_data(3,:,:));
outp_depth1 = squeeze(outp_data(4,:,:));
outp_dtI1 = squeeze(outp_data(5,:,:));

EstPosCoords = sub2ind(size(outp_depth1),EstCoords(1,:),EstCoords(2,:));

EstRads = outp_rad1(EstPosCoords);
EstIs = outp_I1(EstPosCoords);
EstDepths = outp_depth1(EstPosCoords);
EstDtIs = outp_dtI1(EstPosCoords);


TruePars = [EstCoords;EstRads;EstIs;EstCoords;EstDepths;EstDtIs];
TruePars(5,:)=2*TruePars(1,:)*pi()/180;
TruePars(6,:)=2*TruePars(2,:)*pi()/180;
TruePars(1,:)=2*TruePars(1,:)*pi()/180;
TruePars(2,:)=2*TruePars(2,:)*pi()/180;
