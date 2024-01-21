function [alpha0,r0,K0,theta0,phi0] = compute_alldredge_params(original_params)
% Definitions of our original parameters
theta = original_params(1,:);
phi = original_params(2,:);
R = original_params(3,:);
I = original_params(4,:);
r = original_params(7,:);

mu0 = 4*pi*1e-7; % Magnetic permeability of free space
a = 6.378e6; % Earth's radius [m]

% Parameters definitions of Alldredge's method
alpha0 = atan(R./r);
r0 = sqrt(R.^2+r.^2);
K0 = (10^9*mu0*pi*(r0.^2).*(sin(alpha0)).^2.*I)/(a^3);
theta0 = theta;
phi0 = phi;
