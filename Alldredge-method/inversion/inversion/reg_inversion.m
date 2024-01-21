function dp = reg_inversion(looppars, Nmax, fix, loop_nr, ni, subname)
%% Performs a Marquart-Levenberg type non-linear inversion for Alldredge's current loop model
% Input parameters:
%------------------
% looppars: array containing the loop parameters
% Nmax: maximum SH degree used in the forward computations
% fix: set containing the names of loop parameters fixed during the
% estimation e.g. {'K', 'alpha', 'theta', 'phi'};
% loop_nr: number of loop for which the estimation is carried out (0 when
% it is performed on every loop simultaneously)
% ni: total number of iterations passed (for drawing RMS errors)
% subname: name of the function which runs the forward computation (currently: SH_Alldredge)
% Ouptut values:
%------------------
% dp: modifications of the estimated loop parameter

%% HUN
% Az Alldredge-fele hurokparameterek hogyan erkeznek be attol fuggoen, melyeket rogzitjuk 
% Bemeno parameterek:
%--------------------
% looppars: a hurokparameterek tombje Alldredge definicioja szerint
% Nmax: legmagasabb gombi harmonikus fokszam a direktfeladat megoldasokban
% fix: azon hurokparameterek nevei, amelyeket rogzitunk a becsles soran, pl.
% {'K', 'alpha', 'theta', 'phi'}
% loop_nr: annak az aramhuroknak a szama (a 'param' tomb oszlopa), amelyre eppen vegezzuk a becslest (0 when
% ha minden hurokra egyszerre vegezzuk)
% ni: teljes iteracioszam eddig (hibarajzolashoz)
% subname: a fuggveny neve, amely kiszamitja a direktfeladatot (jelenleg: SH_Alldredge)
% Eredmeny:
%--------------------
% dp: a becsulendo parameter modositasa(i) (vektor vagy skalar ha csak egy hurok egy parameterere vegezzuk)
%%

% Initial value of the reguralization factor
% Csillapito tenyezo kezdeti erteke
lam=1e4;

% Maximum possible value of the reguralization constant
% Csillapito tenyezo legnagyobb lehetseges erteke a becsles soran
lammax=1e6;

% Number of estimated loop parameters
% Becsulendo aramhurok parameterek szama
m = length(looppars);

% Using the forward calculation for the data difference vector
% Direktfeladat megoldas, most csak a "mert" es szamitott ertekek kozti
% elteresek megadasahoz
[res, ~]=subname(looppars, Nmax, loop_nr, fix);
res = reshape(res,[],1);
 
% Loading the correct "measured" SHC spectrum
% A helyes ("mert") gombi harmonikus Gauss koefficiensek vektora
ref_spectrum = load('ref_spectrum.mat');
Val = ref_spectrum.ref_spectrum;
Val = reshape(Val,[],1);

% Weights normalizing the data differences and forward computations in sens.m
% Sulyok a "mert"-szamitott elteresek, es a direktfeladat szamitasi eredmenyeinek (sens.m) normalasahoz
w_d = (1e-1+abs(Val));

% Vector of data differences
% "Mert" es szamitott ertekek kozti elteresek vektora
dd=(Val-res);

% Plotting Root Mean Square error
% Legkisebb negyzetes elteres kirajzolasa
RMS=sum(abs(dd))/(length(dd));
figure(1), hold on,
plot(ni, RMS, 'ok')
drawnow

% Factor shrinking parameter modification steps
% Steps are forced to be smaller when smaller RMS values are reached,
% and angular parameters are being estimated
% Szorzo tenyezo a parameter modositasi lepesek lekicsinyitese 
% A szogparameterek becslesenel kisebb lepeseket engedunk csak meg
if (~any(strcmp(fix, 'alpha')))
    %disp('Angle estimation!')
    shrink = 1e-4*RMS;
elseif (~any(strcmp(fix, 'theta'))) || (~any(strcmp(fix, 'phi')))
    %disp('Geogr. coord. estimation!')
    shrink = 1e-4*RMS;
else 
    shrink = 1e-3*RMS;
end

% Computing the Jacobian (sensitivity matrix) S 
% Erzekenyseg (Jacobi) matrix kiszamitasa
S(:,:)=sens(looppars, Nmax, fix, w_d, loop_nr, subname);
dd = dd./w_d;

% Changing the initial regularization factor using the data differences
% Kezdeti csillapito tenyezo megvaltoztatasa az elteresektol fuggoen
lam1=lam*(mean(abs(dd)));
if lam1>lammax
    lam1=lammax;
end
 
% Singular Value Decomposition of the Jacobian S
% Jacobi (erzekenyseg) matrix sajatertek-felbontasa
[U,L,V]=svd(S);
L=diag(L);


% Reconstructing the Jacobian S
% Jacobi (erzekenyseg) matrix felhasznalt SVD komponenseinek megadasa
%--------------------------------------------------------------------
% How much of the singular values we would like to use
% A sajatertekek maximumanak mekkora reszeig szeretnenk oket felhasznalni
SV_bound=1e-4; 
rows = 2*Nmax*Nmax; % total number of rows in S % S sorainak szama

if rows<m
    for nd=1:rows
        if L(nd)>SV_bound*max(L)
            r=nd;
        end
    end
else
    for nd=1:m
        if L(nd)>SV_bound*max(L)
           r=nd;
        end
    end
end
if rows<m
   Ur=U(:,1:rows);
else
   Ur=U(:,1:m);
end
        
Lr=L(1:r);
Vr=V(:,1:r);
Ur=Ur';
%--------------------------------------------------------------------

% Inversion of the corrections
% Parameter modositasok (javitasok) kiszamitasa az SVD inverzzel
if rows<m
   invL2=zeros(r,rows);
else
   invL2=zeros(r,m);
end
    
for ll=1:r                                                                                                                                                    
   invL2(ll,ll)=Lr(ll)/(lam1+Lr(ll)^2);                                                                                                                      
end 


dp=Vr*invL2*Ur*dd;
dp = shrink*dp;
