function [coef, pow] = transform_coefficients_new()
coef = load('coef_new.mat');
pow = load('powers_new.mat');
pow = pow.powers;
coef=coef.coef;
end