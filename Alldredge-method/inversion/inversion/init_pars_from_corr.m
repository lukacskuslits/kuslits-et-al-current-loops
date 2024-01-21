function [alpha0, r0, K0, theta0, phi0] = init_pars_from_corr(perturbations)
correct_pars = load('correct_pars.mat');
correct_looppars = correct_pars.correct_pars;%.looppars;
[alpha0, r0, K0, theta0, phi0] = compute_alldredge_params(correct_looppars);

% add small perturbations

alpha0 = alpha0 + perturbations(1)*alpha0;
r0 = r0 + perturbations(2)*r0;
K0 = K0 + perturbations(3)*K0;
theta0 = theta0 + perturbations(4)*theta0;
phi0 = phi0 + perturbations(5)*phi0;

%Constraining theta between 0.5 and 179 degrees
theta0(theta0>deg2rad(179))=deg2rad(179);
theta0(theta0<deg2rad(0.5))=deg2rad(0.5);

%Constraining alpha between 1 and 75 degrees
alpha0(alpha0>deg2rad(75))=deg2rad(75);
alpha0(alpha0<deg2rad(1))=deg2rad(1);

%Constraining r between Roc and 0.19 Earth radius (a)
a = 6.378e6;
Roc = 3.48e6;
r0(r0>Roc) = Roc;
r0(r0<0.19*a) = a;
end