%% Preparing data
%% Adatok elokeszitese
clear all, close all
Nmax = 15; % Maximum degree of SH expansion
a = 6.378e6; % Earth radius
Roc = 3.48e6; % Outer core radius
nr_of_loops = 5;

%Setting the correct (reference) loop parameters
set_initial_params(nr_of_loops);

% Assign initial parameters using Z field data
% Kezdeti parameterek megadasa kijelolessel 
%[theta0, phi0, alpha0, r0] = load_init_pars();

% Loading the correct loop parameters for computing the reference SH
% spectrum and parameter values
% Helyes hurokparameterek (celertekek) betoltese a referencia (helyes) gombi
% harmonikusokhoz es parameter celertekekhez
pert = [0, 0, 0, 0, 0];
[corr_alpha0, corr_r0, corr_K0, corr_theta0, corr_phi0] = init_pars_from_corr(pert);
%Normalization and logarithms
corr_alpha0 = log(corr_alpha0/(2*pi));
corr_r0 = log(corr_r0/a);
corr_K0 = corr_K0/1e4;
corr_theta0 = log(corr_theta0/(2*pi));
corr_phi0 = log(corr_phi0/(2*pi));


% Assign initial parameters by perturbing the correct values (definitions in alldredge_loop.m)
% Kezdeti parameterek megadasa a helyestol valo elteritessel (definiciojuk az alldredge_loop.m fuggvenyben)
pert = [5e-2, 1e-1, 0, 1e-1, 1e-1];
[alpha0, r0, K0, theta0, phi0] = init_pars_from_corr(pert);

%Normalization and logarithms
alpha0 = log(alpha0/(2*pi));
r0= log(r0/a);
K0 = K0/1e4;
theta0=log(theta0/(2*pi));
phi0=log(phi0/(2*pi));

ini_theta0 = theta0;
ini_phi0 = phi0;
ini_alpha0 = alpha0;
ini_r0 = r0;
ini_K0 = K0;

% Loading the correct original (reference) loop parameters for I 
% Helyes hurokparameterek (celertekek) betoltese - aramerosseg kell innen K
% kiszamitasahoz, mert nem toltjuk meg be
correct_orig_pars_pars = load('correct_pars.mat');
corr_I0 = correct_orig_pars_pars.correct_pars(4,:);

loop_pars_for_ref_spectrum = [corr_alpha0; corr_r0; corr_I0; corr_theta0; corr_phi0]';
fix = 'none';
[ref_spectrum, ref_K] = SH_Alldredge(loop_pars_for_ref_spectrum, Nmax, 0, fix);

% Correct loop parameters for deriving the initial Mean Parameter Error
% (MPE)
% Helyes hurokparameterek a parameterbecsles hibajanak megadasahoz
correct_loop_pars = [corr_alpha0; corr_r0; log(corr_K0); corr_theta0; corr_phi0]';

% Linear inversion of initial K values 
% Kezdeti K normalt magneses momentum megadasa linearis inverzioval
initial_pars = [ini_alpha0; ini_r0; ini_theta0; ini_phi0]';
% A = A_Alldredge(initial_pars, Nmax); % Alakmatrix kiszamitasa
% K0 = alldredge_inv_iter_first(A, ref_spectrum); % LKN inverzio 
% ini_K0 = K0;
% ini_K0= ini_K0'/1e4;

% Calculate initial MPE
% Kezdeti parameterbecslesi hiba megadasa
initial_loop_pars = [ini_alpha0; ini_r0; log(ini_K0); ini_theta0; ini_phi0]';
size_loops = ones(1,nr_of_loops);
%norms = [2*pi*size_loops;a*size_loops;1e4*size_loops;2*pi*size_loops;2*pi*size_loops]';
initial_loop_pars = real(exp(initial_loop_pars));
correct_loop_pars = real(exp(correct_loop_pars));
% ini_MPE = parameter_error_alldredge(correct_loop_pars, initial_loop_pars);
% disp('Initial parameter error:')
% disp(mean(mean(ini_MPE)))

% Saving initial loop parameters and reference SHCs
% Kezdeti hurokparameterek es referencia gombi harmonikus egyutthatok mentese
save('ref_spectrum.mat', 'ref_spectrum')
save('theta.mat', 'theta0')
save('phi.mat', 'phi0')
save('alpha.mat', 'alpha0')
%ref_spectrum = load('ref_spectrum.mat');
%ref_spectrum = ref_spectrum.ref_spectrum;

%loop_pars = [alpha0; r0; correct_looppars(4,:); theta0; phi0]';
%fix = 'none';
%[calc_spectrum, calc_K] = SH_Alldredge(loop_pars, Nmax, fix);


%% Alldredge's inversion iteration
%
figure(1)

% We do not iterate for each loop separately for each loop:
% Ha nem szeretnenk minden forrasra kulon iteraciot:
loop_nr = 0; 

% Start to count total number of iterations for drawing RMS errors
% Iteracios lepesek osszesitett szama az RMS hiba kirajzolasahoz
ni = 1; 

% You can try to do it separately for each loop, the results remain very
% similar:
% Forrasonkent kulon is elvegezheto, az eredmenyek nagyon hasonloak
% maradnak
%for loop_nr = 1:nr_of_loops
for trial = 1:5

% Linear inversion of K values
% K normalt magneses momentumok megadasa linearis inverzioval
A = A_Alldredge(initial_pars, Nmax);
K0 = alldredge_inv_iter_first(A, ref_spectrum);

K0 = K0/1e4; %;
save('K.mat','K0')

est_loop_pars = [alpha0; r0; log(K0'); theta0; phi0];
est_loop_pars = real(exp(est_loop_pars));
MPE1 = parameter_error_alldredge(correct_loop_pars, est_loop_pars');
disp('parameter errors:')
disp(mean(mean(MPE1)))

% Regularised non-linear inversion of r - 3 iterations
% r hurokparameter regularizalt inverzioja - 3 lepes
fix = {'K', 'alpha', 'theta', 'phi'}; %'r_loop', 
subname = @SH_Alldredge;
save('r.mat','r0')
for ii = 1:3
    dp = reg_inversion(r0, Nmax, fix, loop_nr, ni, subname);  
    r0 = real(r0) + real(dp'); %r0.*
    % Constraining r between Roc and 0.19 Earth radius (a)
    % r parameter beszoritasas Rkm es 0.19 Foldsugar (BMH) koze
    r0(a*exp(r0)>Roc) = log(Roc/a);
    r0(a*exp(r0)<0.19*a) = log(0.19);
    ni = ni + 1;
end


save('r.mat','r0') 
initial_pars(:,2)=r0;

% % Regularised non-linear inversion of alpha - 3 iterations
% % alpha hurokparameter regularizalt inverzioja - 3 lepes
% fix = {'K', 'r', 'theta', 'phi'}; %'alpha_loop', 
% subname = @SH_Alldredge;
% save('alpha.mat','alpha0') 
% for ii = 1:3
%     dp = reg_inversion(alpha0, Nmax, fix, loop_nr, ni, subname); 
%     alpha0 = real(alpha0) + real(dp'); % alpha0.*
%     % Constraining alpha between 1 and 75 degrees
%     % alpha hurokparameter beszoritasa 1 es 75 fok koze
%     alpha0(2*pi*exp(alpha0)>deg2rad(75))=log(deg2rad(75)/(2*pi));
%     alpha0(2*pi*exp(alpha0)<deg2rad(1))=log(deg2rad(1)/(2*pi));
%     ni = ni + 1;
% end
% 
% 
% save('alpha.mat','alpha0') 
% initial_pars(:,1)=alpha0;
% 
% theta_orig = theta0;
% % Regularised non-linear inversion of theta0 - 3 iterations
% % theta hurokparameter regularizalt inverzioja - 3 lepes
% fix = {'K', 'alpha', 'r', 'phi'}; %'theta_loop', 
% subname = @SH_Alldredge;
% %save('theta.mat','theta0')
% for ii = 1:3
%     dp = reg_inversion(theta0, Nmax, fix, loop_nr, ni, subname); 
%     theta0 = real(theta0) + real(dp'); %theta0.*
%     % Constraining theta between 0.5 and 179 degrees
%     % theta hurokparameter beszoritasa 0.5 es 179 fok koze
%     theta0(2*pi*exp(theta0)>deg2rad(179))=log(deg2rad(179)/(2*pi));
%     theta0(2*pi*exp(theta0)<deg2rad(0.5))=log(deg2rad(0.5)/(2*pi));
%     ni = ni + 1;
% end
% 
% 
% save('theta.mat','theta0') 
% initial_pars(:,3)=theta0;
% 
% % Regularised non-linear inversion of phi0 - 3 iterations
% % phi hurokparameter regularizalt inverzioja - 3 lepes
% fix = {'K', 'alpha', 'theta', 'r'}; %'phi_loop', 
% subname = @SH_Alldredge;
% %save('phi.mat','phi0')
% for ii = 1:3
%     dp = reg_inversion(phi0, Nmax, fix, loop_nr, ni, subname); 
%     phi0 = real(phi0) + real(dp'); %phi0.*
%     ni = ni + 1;
% end
% 
% save('phi.mat','phi0')
% initial_pars(:,4)=phi0;

% Regularised non-linear inversion of alpha0, theta0, phi0 - 3 iterations
% A szogek regularizalt inverzioja - 3 lepes
fix = {'K', 'r'}; %'phi_loop', 
subname = @SH_Alldredge;
%save('phi.mat','phi0')
loop_pars = zeros(nr_of_loops,3);
for ind_number = 1:nr_of_loops
    loop_pars(ind_number,:) = [alpha0(ind_number), theta0(ind_number), phi0(ind_number)];
end
loop_pars = reshape(loop_pars',[],1);
% loop_pars = [alpha0(1),theta0(1),phi0(1),...
%      alpha0(2),theta0(2),phi0(2),...
%      alpha0(3),theta0(3),phi0(3),...
%      alpha0(4),theta0(4),phi0(4),...
%      alpha0(5),theta0(5),phi0(5)];
for ii = 1:3
    dp = reg_inversion(loop_pars, Nmax, fix, loop_nr, ni, subname); 
    loop_pars = real(loop_pars)+real(dp);
    alpha0 = loop_pars(1:3:3*nr_of_loops);
    alpha0(2*pi*exp(alpha0)>deg2rad(75))=log(deg2rad(75)/(2*pi));
    alpha0(2*pi*exp(alpha0)<deg2rad(1))=log(deg2rad(1)/(2*pi));
    loop_pars(1:3:3*nr_of_loops) = alpha0;
    theta0 = loop_pars(2:3:3*nr_of_loops);
    theta0(2*pi*exp(theta0)>deg2rad(179))=log(deg2rad(179)/(2*pi));
    theta0(2*pi*exp(theta0)<deg2rad(0.5))=log(deg2rad(0.5)/(2*pi));
    loop_pars(2:3:3*nr_of_loops) = theta0;
    ni = ni + 1;
end

alpha0 = loop_pars(1:3:3*nr_of_loops)';
theta0 = loop_pars(2:3:3*nr_of_loops)';
phi0 = loop_pars(3:3:3*nr_of_loops)';
save('alpha.mat','alpha0') 
initial_pars(:,1)=alpha0;
save('theta.mat','theta0')
initial_pars(:,3)=theta0;
save('phi.mat','phi0')
initial_pars(:,4)=phi0;

est_loop_pars = [alpha0; r0; log(K0'); theta0; phi0];
est_loop_pars = real(exp(est_loop_pars));
MPE2 = parameter_error_alldredge(correct_loop_pars, est_loop_pars');
disp('parameter errors:')
disp(mean(mean(MPE2)))

end
%end

% Calculating final MPE
% Vegso parameterbecslesi hiba megadasa
est_loop_pars = [alpha0; r0; log(K0'); theta0; phi0];
est_loop_pars = real(exp(est_loop_pars));
MPE3 = parameter_error_alldredge(correct_loop_pars, est_loop_pars');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial parameter estimation error: 5 loops
% Provided that \delta (alpha0,K0) = 5%; \delta (r0) = 10%
% theta0: 1e-2; phi0: 1e-2; ini_MPE: 0.0697 --> MPE: 0.0689 | decr_ratio = -0.115 
% theta0: 5e-2; phi0: 5e-2; ini_MPE: 0.1272 --> MPE: 0.1126 | decr_ratio = -0.115
% theta0: 1e-1; phi0: 1e-1; ini_MPE: 0.2211 --> MPE: 0.2257 | decr_ratio = 0.0208 

% Initial parameter estimation error: 10 loops
% Provided that \delta (alpha0,K0) = 5%; \delta (r0) = 10%
% theta0: 5e-2; phi0: 5e-2; ini_MPE: 0.2691 --> MPE: 0.1797 | decr_ratio = -0.332
% theta0: 1e-1; phi0: 1e-1; ini_MPE: 0.4079 --> MPE: 0.3589 | decr_ratio = -0.12 

% Initial parameter estimation error: 20 loops
% Provided that \delta (alpha0,K0) = 5%; \delta (r0) = 10%
% theta0: 5e-2; phi0: 5e-2; ini_MPE: 0.3475 --> MPE: 0.4360  | decr_ratio = 0.255
% theta0: 1e-1; phi0: 1e-1; ini_MPE: 1.0628 --> MPE: 1.0837 | decr_ratio = 0.0197 