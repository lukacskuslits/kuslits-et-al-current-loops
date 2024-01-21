function val = gaussC(max_val, x, y, sigma, center)
%% Calculates the values of a 2d Gaussian distribution mask around a coordinate center
%Input params:
%-------------
% max_val: value at the center (1 when data are normalized)
% x: input x coordinates
% y: input y coordinates
% sigma: standard deviation
% center: coordinate values of the center of the distribution (xc, yc)
%Ouptut values:
%--------------
% val: vector containing gaussian distribution values

xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
val       = max_val*exp(-exponent); 
