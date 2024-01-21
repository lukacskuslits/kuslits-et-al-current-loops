function [dgnm_dalpha_anal, dgnm_dalpha_num] = test_sens(n,m)
alpha0 = pi/50; 
r0 = 3.1e6; 
I0 = 1e7; 
theta0 = pi/3; 
phi0 = pi/3;
corr_loop_pars = [alpha0; r0; I0; theta0; phi0]';
fix = 'none';
Nmax = 15;
subname = @SH_Alldredge;

[~, K0] = SH_Alldredge(corr_loop_pars, Nmax, 0, fix);

save('r.mat', 'r0')
save('K.mat', 'K0')
save('theta.mat', 'theta0')
save('phi.mat', 'phi0')
save('alpha.mat', 'alpha0')

dgnm_dalpha_anal = analytic_sol(n,m,alpha0,r0,K0,theta0,phi0);
dgnm_dalpha_num = numeric_sol(alpha0, subname);

end


function dgnm_dalpha = analytic_sol(n,m,alpha,r,K,theta,phi)

a = 6.378e6; %Earth's radius [m]

numerator_11 = sqrt(2)*K*(r/a)^(n-1);
numerator_12 = n*cos(alpha)*evaluate_legendre(n,1,alpha)...
               -(n+1)*evaluate_legendre(n-1,1,alpha);
numerator_13 = cos(m*phi)*evaluate_legendre(n,m,theta);
denominator_1 = sqrt(n*(n+1))*((cos(alpha))^2-1);

numerator_21 = numerator_11;
numerator_22 = cos(alpha)*cos(m*phi)*evaluate_legendre(n,1,alpha);
numerator_23 = evaluate_legendre(n,m,theta);
denominator_2 = sqrt(n*(n+1))*(sin(alpha))^2;

dgnm_dalpha_1 = (numerator_11*numerator_12*numerator_13)/denominator_1;
dgnm_dalpha_2 = (numerator_21*numerator_22*numerator_23)/denominator_2;

dgnm_dalpha = -dgnm_dalpha_1 - dgnm_dalpha_2;
end


function dgnm_dalpha = numeric_sol(alpha, subname)

fix = {'r', 'K', 'theta','phi'};
pert=1e-4; %Small perturbation

dp = abs(alpha*pert);
Nmax = 15;

%Preallocating Jacobian
%S=zeros(rows,1);

%w_p = looppar;
  
%Loop calculating S element-wise
pnew=alpha+dp; 
pnew2=alpha-dp;

% %Sliced parameters for the forward calc. functions
% %Calling forwards using the unchanged and changed params - a
[Fmvj, ~]=subname(pnew,Nmax,0,fix);
Fmvj=reshape(Fmvj,[],1);
[Fmvvj, ~]=subname(pnew2,Nmax,0,fix);
Fmvvj=reshape(Fmvvj,[],1);

%Calculating S columnwise
delta_d = (Fmvj-Fmvvj);
dgnm_dalpha=delta_d/(2*dp);

end