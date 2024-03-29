function induced = ind_approx_new(fi, r, depth)
% fi = 0;%(pi/3)*(pi^(-1));
% depth = 1.25e-1;%2.5e5/8e5;
% r = 1.54e-1;%1e5/6.5e5;
coef = load('coef_use.mat');
coef = coef.coef;
pow = load('powers_use.mat');
pow = pow.powers;
X = [fi, r, depth];
XX = repmat(X,[1,1,length(pow(1,1,:))]);
X_pow = XX.^pow;%bsxfun(@power,XX,pow);
poly_x = prod(X_pow,2);
poly_x = squeeze(poly_x);
induced = 2.5e-3 - exp(-5.095 + poly_x * coef');
end