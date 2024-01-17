function induced = ind_approx_orig(fi, depth, r)
% fi = 0;%(pi/3)*(pi^(-1));
% depth = 1.25e-1;%2.5e5/8e5;
% r = 1.54e-1;%1e5/6.5e5;
X = [fi, depth, r];
coef = load('coef.mat');
pow = load('powers.mat');
pow = pow.powers;
coef=coef.coef;
pow = cast(pow, 'like', coef);
XX = repmat(X,[1,1,length(pow)]);
pow = repmat(pow,[1,1, length(XX)]);
pow = permute(pow, [3,2,1]);
% disp(size(X))
X_pow = XX.^pow;%bsxfun(@power,XX,pow);
poly_x = prod(X_pow,2);
poly_x = squeeze(poly_x);
induced = poly_x * coef';
end
