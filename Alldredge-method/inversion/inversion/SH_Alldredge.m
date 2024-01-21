function [spectrum, Ks] = SH_Alldredge(loop_pars, Nmax, loop_nr, fix)
%% Calculates the spherical harmonic coefficients (SHC) of the magnetic
%% field of loop_nr current loops.
% Input parameters:
%------------------
% loop_pars: array containing Alldredge's version of loop parameters
% Nmax: maximum degree of SH expansion
% loop_nr: number of the loop (column of array 'param') for which the estimation is carried out (0 when
% it is performed on every loop simultaneously)
% fix: set containing the names of loop parameters fixed during the
% estimation e.g. {'K', 'alpha', 'theta', 'phi'}
% Output values:
%------------------
% spectrum: array containing SHCs
% Ks: vector of the normalised magnetic moments of each loop

%% HUN
% Kiszamitja a Gauss koefficienseket loop_nr aramhurok eseten
% Bemeno parameterek:
%--------------------
% loop_pars: a hurokparameterek tombje Alldredge definicioja szerint
% Nmax: legmagasabb gombi harmonikus fokszam
% loop_nr: annak az aramhuroknak a szama (a 'param' tomb oszlopa), amelyre eppen vegezzuk a becslest (0 when
% ha minden hurokra egyszerre vegezzuk)
% fix: azon hurokparameterek nevei, amelyeket rogzitunk a becsles soran, pl.
% {'K', 'alpha', 'theta', 'phi'}
% Eredmeny:
%--------------------
% spectrum: Gauss egyutthatok tombje
% Ks: valamennyi hurok normalt magneses momentumat tartalmazo vektor
%%

sloops = size(loop_pars);

% fixed parameter values are red from the corresponding .mat
% files
if strcmp(fix,'none')
    nloops = sloops(1);
else
    nloops = sloops(1)*sloops(2);
end
      
if any(strcmp(fix, 'alpha'))
   alpha = load('alpha.mat');
   alpha0 = alpha.alpha0;
end
if any(strcmp(fix, 'K'))
    K = load('K.mat');
    K0 = K.K0;
end
if any(strcmp(fix, 'r'))
    r = load('r.mat');
    r0 = r.r0;
end
if any(strcmp(fix, 'theta'))
    theta = load('theta.mat');
    theta0 = theta.theta0;
end
if any(strcmp(fix, 'phi'))
    phi = load('phi.mat');
    phi0 = phi.phi0;
end


spectrum = zeros(Nmax,Nmax,2);
Ks = zeros(1,nloops);
% Vakuum magneses permeabilitasa
mu0 = 4*pi*1e-7; % Magnetic permeability of free space [H/m] 
% Foldsugar
a = 6.378e6; % Earth's radius [m]

%try to invert the angles in the same step
if length(fix)<4
    niter = nloops/3;
else
    niter = nloops;
end

for ii = 1:niter
    if ~any(strcmp(fix, 'alpha'))
         par_nr = 1;
         if length(fix)<4
             element_nr = ii*3-2;
         else
             element_nr = ii;
         end
         alpha = structure_input('alpha', loop_pars, element_nr, loop_nr, par_nr, nloops, fix);
    else
        alpha=alpha0(ii);
    end
    alpha = 2*pi*exp(alpha);
    %
    if ~any(strcmp(fix, 'r'))
        par_nr = 2;
        r = structure_input('r', loop_pars, ii, loop_nr, par_nr, nloops, fix);
%     elseif length(fix)<4
%         [~,col] = ind2sub([5,5],ii);
%         r=r0(col);
    else
        r=r0(ii);
    end
    r= a*exp(r);
    %
    if ~any(strcmp(fix, 'K'))
        I=loop_pars(ii,3); 
        K = (10^9*mu0*pi*r^2*(sin(alpha))^2*I)/(a^3);
    elseif strcmp(fix, 'K_loop')
        K0 = load('K.mat');
        K0 = K0.K0;
        if loop_nr == 1
           K = [loop_pars(ii), K0(loop_nr+1:nloop)];
        elseif loop_nr == nloop
           K=[K0(1:loop_nr), loop_pars(ii)];
        else
           K=[K0(1:loop_nr), loop_pars(ii), K0(loop_nr+1:nloop)];
        end
%     elseif length(fix)<4
%         [~,col] = ind2sub([5,5],ii);
%         K=K0(col);
%         K= 10^4*K;
    else
        K = K0(ii);
        K= 10^4*K;
    end
    %
    if ~any(strcmp(fix, 'theta'))
        par_nr = 4;
        if length(fix)<4
            element_nr = ii*3-1;
        else
            element_nr = ii;
        end
        theta = structure_input('theta', loop_pars, element_nr, loop_nr, par_nr, nloops, fix);
    else 
        theta = theta0(ii);
    end
    %
    if ~any(strcmp(fix, 'phi'))
        par_nr = 5;
        if length(fix)<4
            element_nr = ii*3;
        else
            element_nr = ii;
        end
        phi = structure_input('phi', loop_pars, element_nr, loop_nr, par_nr, nloops, fix);
    else
        phi = phi0(ii);
    end
    %
    %Normalization and exponentials
    theta=2*pi*exp(theta);
    phi=2*pi*exp(phi);
    
    %Saving Ks
    Ks(ii) = K;
    for degree=1:Nmax
        for order=1:degree
            [gnm, hnm] = alldredge_loop(degree, order, K, alpha, r, theta, phi);
            spectrum(degree, order, 1) = spectrum(degree, order, 1)+gnm;
            spectrum(degree, order, 2) = spectrum(degree, order, 2)+hnm;
        end
    end
end
end