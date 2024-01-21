function induced = ind_approx_single(fi, r, depth)
%% Computes an approximate induced contribution above a single radially aligned loop
%Input paramas:
%--------------
%fi: matrix containing angular distances between the points of observation and the loop centre
%r: matrix containing the normalized loop radius
%depth: matrix containing the normalized depth of the loop calculated from the CMB
%Input files:
%------------
% coef_new.mat: matrix containing polynomial coefficients used for the
    % approximations of diffusion screening
% powers_new.mat: matrix containing polynomial powers used for the
    % approximations of diffusion screening
%Output values:
%--------------
% induced: induced radial magnetic field value [T]

X = [fi, r, depth];

%------ induction single fi------------------------------------
coef = load('coef_new.mat');
pow = load('powers_new.mat');
coef = coef.coef;
pow = pow.powers;
pow = cast(pow, 'double');
XX = repmat(X,[1,1,length(pow)]);
XX = squeeze(XX);
XX = permute(XX,[2,1]);
X_pow = XX.^pow;
poly_x = prod(X_pow,2);
induced = 2.5e-3 - exp(-5.095 + coef*poly_x);

