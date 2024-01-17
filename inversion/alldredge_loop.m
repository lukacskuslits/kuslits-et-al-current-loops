function [gnm, hnm] = alldredge_loop(n, m, K, alpha, r, theta_0, phi_0)
%% Computes the sine and cosine SHCs of order n, degree m for a single current loop
% Input parameters:
%-----------------------
% n: SH degree
% m: SH order
% K: normalized magnetic moment
% alpha: half angle of sight from EC to the loop
% r: distance from EC to the edge of the loop
% theta_0: spherical colatitude
% phi_0: spherical longitude
% Output values:
%-----------------------
% gnm sine SHC
% hnm cosine SHC

%%HUN
% Szinusz es koszinusz gombi harmonikusok Gauss koefficienseit adja meg egy darab aramhurok eseten
% Bemeno parameterek:
%--------------------
% n: fok
% m: rend
% K: normalt magneses momentum
% alpha: a hurok Foldkozeppontbol vett latoszogenek fele [radian]
% r: Foldkozeppont es a hurokaram szele kozti tavolsag [m]
% theta_0: gombi segedszelesseg [radian]
% phi_0: gombi hosszusag [radian]
% Eredmeny:
%--------------------
% gnm szinusz egyutthatok
% hnm koszinusz egyutthatok
%%


a = 6.378e6; % Earth's radius [m]

x = cos(alpha);
legendre_1 = legendre(n, x);
legendre_1 = legendre_1(2);
g0n = K*((2*n)/(sin(alpha)))*((r/a)^(n-1))*(legendre_1/sqrt(2*n*(n+1)));

x_m = cos(theta_0);
legendre_m = legendre(n, x_m, 'sch');
legendre_m = legendre_m(m);
gnm = g0n*legendre_m*cos(m*phi_0);
hnm = g0n*legendre_m*sin(m*phi_0);


end