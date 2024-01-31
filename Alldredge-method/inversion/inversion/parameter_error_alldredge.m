function RPE = parameter_error_alldredge(true_loops, est_loops)
%% Calculates the relative parameter errors for a given estimation of the loop parameters 
%% using the Alldredge method
%%(Finding the physically nearest true loop is not needed here as initial
%%loops were fixed to be the same during the stability analysis)
% Input params:
%--------------
%est_loops: array containing the estimated loop parameters
%true_loops: array containing the true loop parameters
% Output values:
%---------------
%RPE: array containing normalized relative parameter estimation errors

[scale_alpha, scale_r, scale_K, scale_theta, scale_phi, Kmax] = get_scaling_values();
scaling_vals = [scale_alpha; scale_r; scale_K; scale_theta; scale_phi];
true_loops = bsxfun(@times,scaling_vals, true_loops);
est_loops = bsxfun(@times,scaling_vals, est_loops);

%Tolerance to avoid very large relative errors in case of small normalized magnetic moment values (the magnetic moment values can cross zero)
tolerance = [0;0;Kmax/100;0;0]; 

RPE = abs(true_loops-est_loops)./bsxfun(@plus,tolerance, abs(true_loops));

end

function [scale_alpha, scale_r, scale_K, scale_theta, scale_phi, Kmax] = get_scaling_values()

Re = 6.378e6;
mu0 = 4*pi*1e-7;
rmax = 9e5;
Imax = 1e9;
alpha_max = deg2rad(75);

%maximum possible value of the normalized magnetic moment as a base of tolerance
Kmax = (10^9*mu0*pi*rmax^2*(sin(alpha_max))^2*Imax)/(Re^3);
scale_alpha = 2*pi;
scale_theta = 2*pi;
scale_phi = 2*pi;
scale_r = Re;
scale_K = 1e4;

end
