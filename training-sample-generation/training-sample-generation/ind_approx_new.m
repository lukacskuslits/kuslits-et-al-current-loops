function induced = ind_approx_new(fi, r, depth)
%% Computes an approximate induced contribution for a radially aligned loop
%Input paramas:
%--------------
%fi: matrix containing angular distances between the points of observation and the loop centre
%depth: matrix containing the normalized depth of the loop calculated from the CMB
%a: matrix containing the normalized loop radius
%Input files:
%------------
% coef_use.mat: preprocessed (transformed) matrix containing polynomial coeffitients used for the
    % approximations of diffusion screening
% powers_use.mat: preprocessed (transformed) matrix containing polynomial powers used for the
    % approximations of diffusion screening
%Output values:
%--------------
% induced: vector of induced radial magnetic field data [T]

coef = load('coef_use.mat');
coef = coef.coef;
pow = load('powers_use.mat');
pow = pow.powers;
X = [fi, r, depth];
XX = repmat(X,[1,1,length(pow(1,1,:))]);
X_pow = XX.^pow;
poly_x = prod(X_pow,2);
poly_x = squeeze(poly_x);
induced = 2.5e-3 - exp(-5.095 + poly_x * coef');
end
