function transform_coefficients_new(deg_res)
%% Transforms the matrices of exponential function-based polynomial coefficients to be used for estimating diffusion
% Input params:
%--------------
% deg_res: resolution of output field maps in spherical coordinates (degrees)
%Input files:
%------------
% coef_new.mat: matrix containing polynomial coefficients used for the
    % approximations of diffusion screening
% powers_new.mat: matrix containing polynomial powers used for the
    % approximations of diffusion screening
%Output files:
%-------------
% coef_use.mat: preprocessed (transformed) matrix containing polynomial coefficients used for the
    % approximations of diffusion screening
% powers_use.mat: preprocessed (transformed) matrix containing polynomial powers used for the
    % approximations of diffusion screening

coef = load('coef_new.mat');
coef = coef.coef;
powers = load('powers_new.mat');
powers = powers.powers;
lat_dim = 180/deg_res;
long_dim = 360/deg_res;

X = ones(lat_dim*long_dim,3);
powers = cast(powers, 'double');
XX = repmat(X,[1,1,length(powers)]);
powers = repmat(powers,[1,1, length(XX)]);
powers = permute(powers, [3,2,1]);
save('coef_use.mat', 'coef')
save('powers_use.mat', 'powers')

end
