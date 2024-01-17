function [coef, powers] = transform_coefficients_new(lat_dim, long_dim)
% dt = 25*3.15e7;
% mu = 4*pi*1e-7;
coef = load('coef_new.mat');
coef = coef.coef;
powers = load('powers_new.mat');
powers = powers.powers;
%lat_dim = 96;%180/deg_res;
%long_dim = 192;%360/deg_res;

X = ones(lat_dim*long_dim,3);
%powers = cast(powers, 'like', coef);
powers = cast(powers, 'double');
XX = repmat(X,[1,1,length(powers)]);
powers = repmat(powers,[1,1, length(XX)]);
powers = permute(powers, [3,2,1]);
save('coef_use.mat', 'coef')
save('powers_use.mat', 'powers')
% pow = pow.powers;
% coef=coef.coef;
end