function [coef, pow] = transform_coefficientsPars()
coef = load('coef.mat');
pow = load('powers.mat');
pow = pow.powers;
coef=coef.coef;
X = ones(4050,3);
pow = cast(pow, 'like', coef);
XX = repmat(X,[1,1,length(pow)]);
pow = repmat(pow,[1,1, length(XX)]);
pow = permute(pow, [3,2,1]);
end