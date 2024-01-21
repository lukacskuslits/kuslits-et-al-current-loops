function A = A_Alldredge(loop_pars, Nmax)
%% Computes the systems matrix A for K
%% K-hoz kiszamitja a linearis egyenletrendszer alakmatrixat

sloops = size(loop_pars);
nloops = sloops(1);

spectrum = zeros(Nmax,Nmax,2);
%mu0 = 4*pi*1e-7; % Magnetic permeability of free space [H/m]
a = 6.378e6; % Earth's radius [m]

A = zeros(Nmax*Nmax*2,nloops);
for ii = 1:nloops
    alpha=loop_pars(ii,1);
    r=loop_pars(ii,2);
    theta_0=loop_pars(ii,3);
    phi_0=loop_pars(ii,4);
    
    %Normalization and exponentials
    alpha = 2*pi*exp(alpha);
    r= a*exp(r);
    theta_0=2*pi*exp(theta_0);
    phi_0=2*pi*exp(phi_0);

    K = 1; %(10^9*mu0*pi*r^2*(sin(alpha))^2*I)/(a^3);
    for degree=1:Nmax
        for order=1:degree
            [gnm, hnm] = alldredge_loop(degree, order, K, alpha, r, theta_0, phi_0);
            spectrum(degree, order, 1) = gnm;
            spectrum(degree, order, 2) = hnm;
        end
    end
    A(:,ii) = reshape(spectrum,[],1);
end
end