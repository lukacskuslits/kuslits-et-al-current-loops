function induced = ind_approx_single(fi, r, depth)
X = [fi, r, depth];

%------ induction single fi------------------------------------
coef = load('coef_new.mat');
pow = load('powers_new.mat');
coef = coef.coef;
pow = pow.powers;
%pow = cast(pow, 'like', coef);
pow = cast(pow, 'double');
XX = repmat(X,[1,1,length(pow)]);
XX = squeeze(XX);
XX = permute(XX,[2,1]);
%powr = repmat(pow,[1,1, length(XX)]);
%powr = permute(powr, [3,2,1]);
X_pow = XX.^pow;
poly_x = prod(X_pow,2);
%poly_x = squeeze(poly_x);
induced = 2.5e-3 - exp(-5.095 + coef*poly_x);
%induced = coef * poly_x;