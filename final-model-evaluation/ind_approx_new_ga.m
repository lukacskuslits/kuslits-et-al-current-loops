function induced = ind_approx_new_ga(fi, r, depth, coef, pow)
% fi = 0;%(pi/3)*(pi^(-1));
% depth = 1.25e-1;%2.5e5/8e5;
% r = 1.54e-1;%1e5/6.5e5;
X = [fi, r, depth];
XX = repmat(X,[1,1,length(pow(1,1,:))]);
X_pow = XX.^pow;%bsxfun(@power,XX,pow);
poly_x = prod(X_pow,2);
poly_x = squeeze(poly_x);
induced = 2.5e-3 - exp(-5.095 + poly_x * coef');
end