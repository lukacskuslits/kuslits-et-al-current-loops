%% Preparing data
%% Adatok elokeszitese
clear all, close all
Nmax = 15; % Maximum degree of SH expansion
a = 6.378e6; % Earth radius
Roc = 3.48e6; % Outer core radius

% Loading the correct (reference) loop parameters
% Helyes hurokparameterek (celertekek) betoltese
correct_pars = load('correct_pars.mat');
correct_looppars = correct_pars.correct_pars.looppars;
correct_looppars(1,:)=rad2deg(correct_looppars(1,:));
correct_looppars(2,:)=rad2deg(correct_looppars(2,:));
correct_Ks = correct_pars.correct_pars.Ks;

% Assign initial parameters using Z field data
% Kezdeti parameterek megadasa kijelolessel 
%[theta0, phi0, alpha0, r0] = load_init_pars();

% Assign initial parameters by perturbing the correct values (definitions in alldredge_loop.m)
% Kezdeti parameterek megadasa a helyestol valo elteritessel (definiciojuk az alldredge_loop.m fuggvenyben)
pert = [1e-3,1e-3,0,0];
[theta0, phi0, alpha0, r0] = init_pars_from_corr(pert);
ini_theta0 = theta0;
ini_phi0 = phi0;
ini_alpha0 = alpha0;
ini_r0 = r0;

% Loading the correct loop parameters for computing the reference SH spectrum
% Helyes hurokparameterek (celertekek) betoltese a referencia (helyes) gombi
% harmonikusokhoz
pert = [0;0;0;0];
[corr_theta0, corr_phi0, corr_alpha0, corr_r0] = init_pars_from_corr(pert);
corr_loop_pars = [corr_alpha0; corr_r0; correct_looppars(4,:); corr_theta0; corr_phi0]';
fix = 'none';
[ref_spectrum, ref_K] = SH_Alldredge(corr_loop_pars, Nmax, 0, fix);

% Correct loop parameters for deriving the initial Mean Parameter Error
% (MPE)
% Helyes hurokparameterek a parameterbecsles hibajanak megadasahoz
correct_loop_pars = [corr_alpha0; corr_r0; correct_Ks; corr_theta0; corr_phi0]';

% Linear inversion of initial K values 
% Kezdeti K normalt magneses momentum megadasa linearis inverzioval
initial_pars = [alpha0; r0; theta0; phi0]';
A = A_Alldredge(initial_pars, Nmax); % Alakmatrix kiszamitasa
K0 = alldredge_inv_iter_first(A, ref_spectrum); % LKN inverzio 
ini_K0 = K0;

% Calculate initial MPE
% Kezdeti parameterbecslesi hiba megadasa
initial_loop_pars = [ini_alpha0; ini_r0; ini_K0'; ini_theta0; ini_phi0]';
ini_MPE = parameter_error_alldredge(correct_loop_pars, initial_loop_pars);

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
%for loop_nr = 1:5
for trial = 1:10

% Linear inversion of K values
% K normalt magneses momentumok megadasa linearis inverzioval
A = A_Alldredge(initial_pars, Nmax);
K0 = alldredge_inv_iter_first(A, ref_spectrum);

%K0 = calc_K;
ini_K0 = K0;
save('K.mat','K0')

% Regularised non-linear inversion of r - 3 iterations
% r hurokparameter regularizalt inverzioja - 3 lepes
fix = {'K', 'alpha', 'theta', 'phi'}; %'r_loop', 
subname = @SH_Alldredge;
save('r.mat','r0')
for ii = 1:3
    dp = reg_inversion(r0, Nmax, fix, loop_nr, ni, subname);  
    r0 = r0 - r0.*dp';
    % Constraining r between Roc and 0.19 Earth radius (a)
    % r parameter beszoritasas Rkm es 0.19 Foldsugar (BMH) koze
    r0(r0>Roc) = Roc;
    r0(r0<0.19*a) = a;
    ni = ni + 1;
end


save('r.mat','r0') 
initial_pars(:,2)=r0;

% Regularised non-linear inversion of alpha - 3 iterations
% alpha hurokparameter regularizalt inverzioja - 3 lepes
fix = {'K', 'r', 'theta', 'phi'}; %'alpha_loop', 
subname = @SH_Alldredge;
save('alpha.mat','alpha0') 
for ii = 1:3
    dp = reg_inversion(alpha0, Nmax, fix, loop_nr, ni, subname); 
    alpha0 = alpha0 - alpha0.*dp';
    % Constraining alpha between 1 and 75 degrees
    % alpha hurokparameter beszoritasa 1 es 75 fok koze
    alpha0(alpha0>deg2rad(75))=deg2rad(75);
    alpha0(alpha0<deg2rad(1))=deg2rad(1);
    ni = ni + 1;
end


save('alpha.mat','alpha0') 
initial_pars(:,1)=alpha0;

theta_orig = theta0;
% Regularised non-linear inversion of theta0 - 3 iterations
% theta hurokparameter regularizalt inverzioja - 3 lepes
fix = {'K', 'alpha', 'r', 'phi'}; %'theta_loop', 
subname = @SH_Alldredge;
%save('theta.mat','theta0')
for ii = 1:3
    dp = reg_inversion(theta0, Nmax, fix, loop_nr, ni, subname); 
    theta0 = theta0 - theta0.*dp';
    % Constraining theta between 0.5 and 179 degrees
    % theta hurokparameter beszoritasa 0.5 es 179 fok koze
    theta0(theta0>deg2rad(179))=deg2rad(179);
    theta0(theta0<deg2rad(0.5))=deg2rad(0.5);
    ni = ni + 1;
end


save('theta.mat','theta0') 
initial_pars(:,3)=theta0;

% Regularised non-linear inversion of phi0 - 3 iterations
% phi hurokparameter regularizalt inverzioja - 3 lepes
fix = {'K', 'alpha', 'theta', 'r'}; %'phi_loop', 
subname = @SH_Alldredge;
%save('phi.mat','phi0')
for ii = 1:3
    dp = reg_inversion(phi0, Nmax, fix, loop_nr, ni, subname); 
    phi0 = phi0 - phi0.*dp';
    ni = ni + 1;
end

save('phi.mat','phi0')
initial_pars(:,4)=phi0;
end
%end

% Calculating final MPE
% Vegso parameterbecslesi hiba megadasa
est_loop_pars = [alpha0; r0; K0'; theta0; phi0]';
MPE = parameter_error_alldredge(correct_loop_pars, est_loop_pars);

%%Systhematic test results:
% Initial parameter estimation error: 
% theta0: 1e-3; alpha0: 1e-2; ini_MPE: 0.0028 --> MPE: 0.0027
% theta0: 1e-2; alpha0: 1e-2; ini_MPE: 0.0045 --> MPE: 0.0071
% theta0: 1e-3; r0: 1e-2; ini_MPE: 0.0104 --> MPE: 0.0089
% theta0: 1e-2; r0: 1e-2; ini_MPE: 0.0130 --> MPE: 0.0088
% theta0: 1e-3; phi0: 1e-3; ini_MPE: 0.0011 --> MPE: 0.0013
% theta0: 1e-3; phi0: 1e-2; ini_MPE: 0.0096 --> MPE: 0.02
% theta0: 1e-2; phi0: 1e-2; ini_MPE: 0.0114 --> MPE: 0.0230
