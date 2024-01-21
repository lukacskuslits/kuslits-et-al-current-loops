function BRR = total_field_approx_new(pars, long_size, lat_size)
%% Generates the total radial magnetic field of a loop model
%Input params:
%-------------
% pars: matrix containing all the parameters of every loop in the model
% long_size: number of longitude bins to compute the radial MF for
% lat_size: number of latitude bins to compute the radial MF for 
%Output values:
%--------------
%BRR: vector containing radial magnetic field data [T]

Roc = 3.478e6;
mu = 4 * pi * 1e-7;
theta=0:pi/(long_size-1):pi();
lambda=0:2*pi/(lat_size-1):2*pi();


[lambda,theta]=meshgrid(lambda,theta);
lambda = lambda';
theta = theta';
lambda=reshape(lambda,numel(lambda),1);
theta=reshape(theta,numel(theta),1);
BRR = zeros(length(theta),1);

for kl=1:length(pars(1,:))
thetl = pars(1,kl);
laml = pars(2,kl);
a = pars(3,kl);
I = pars(4,kl);
thet = pars(5,kl); 
lam = pars(6,kl);
rloop = pars(7,kl);
dtI = pars(8,kl);

dx = rloop * sin(thetl) * cos(laml);
dy = rloop * sin(thetl) * sin(laml);
dz = rloop * cos(thetl);
% 

Rocx = Roc * sin(theta) .* cos(lambda);
Rocy = Roc * sin(theta) .* sin(lambda);
Rocz = Roc * cos(theta);

% BR at CMB
x = Rocx - dx;
y = Rocy - dy;
z = Rocz - dz;

%%
% perform coordinate (position) rotation
t = thet;
l = lam;
alf = 0;
R1=[1 0 0; 0 cos(alf) -sin(alf); 0 sin(alf) cos(alf)]; 
R2=[cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
R3=[cos(l) -sin(l) 0; sin(l) cos(l) 0; 0 0 1];
rot = R3*R2*R1;
roti = inv(rot);

BB = [x'; y'; z'];
LX = length(x);
BBr = zeros(3,LX,3);

for ii = 1:3
BBr(ii,:,1) = BB(1,:)*roti(ii);
BBr(ii,:,2) = BB(2,:)*roti(ii+3);
BBr(ii,:,3) = BB(3,:)*roti(ii+6);

end
BBr = sum(BBr,3);


x = BBr(1,:)';
y = BBr(2,:)';
z = BBr(3,:)';
%%

r = sqrt(x.^2 + y.^2 + z.^2);
rho = sqrt(x.^2 + y.^2);

[Bx, By, Bz] = loop_analytic(x, y, z, a, r, rho, I, mu);


%%
% perform rotation of vector field Descartes components
BB = [Bx'; By'; Bz'];
BBr = zeros(3,LX,3);
%jj = 1;
for ii = 1:3
BBr(ii,:,1) = BB(1,:)*rot(ii);
BBr(ii,:,2) = BB(2,:)*rot(ii+3);
BBr(ii,:,3) = BB(3,:)*rot(ii+6);

end
BBr = sum(BBr,3);

Bx = BBr(1,:)';
By = BBr(2,:)';
Bz = BBr(3,:)';


dBRR1 = Bx .* sin(theta) .* cos(lambda);
dBRR2 = By .* sin(theta) .* sin(lambda);
dBRR3 = Bz .* cos(theta);
dBRR = dBRR1 + dBRR2 + dBRR3;
%% Compute the induced field components if there is a time variation in the loop current
if dtI~=0
    lRoc = Rocx.^2 + Rocy.^2 + Rocz.^2;
    lloop = dx.^2 + dy.^2 + dz.^2;
    fi = acos((Rocx.*dx + Rocy.*dy + Rocz.*dz)./sqrt(lRoc.*lloop))/pi;
    depth = (Roc-rloop)/9e5;
    a = a/8e5;%6.5e5;
    depth = repmat(depth, length(fi), 1);
    a = repmat(a, length(fi), 1);
    dBRR_i = ind_approx_new(fi, a, depth); 
    dBRR = dBRR + dBRR_i*dtI;
end
BRR = BRR + dBRR;
end
end


function [Bx, By, Bz] = loop_analytic(x, y, z, a, r, rho, I, mu) 

%% Computes the cartesian magnetic field components of a single current loop
%Input params:
%x: cartesian x distances between the points of observation and the centre
%of the loop
%y: cartesian y distances between the points of observation and the centre
%of the loop
%z: cartesian z distances between the points of observation and the centre
%of the loop
%a: loop radius
%r: euclidean distances between the points of observation and the centre
%of the loop
%rho: cylindrical axial distances between the points of observation and the centre
%of the loop

%Expressions, after Simpson et al. (2001):
% alpha^2 = a^2 + r^2 - 2 a rho
% beta^2 = a^2 + r^2 + 2 a rho
% k^2 = 1 - alpha^2 beta^(-2)  
% C = mu I pi^(-1)

alpha_2 = a^2 + r.^2 - 2 * a * rho;  
beta_2 = a^2 + r.^2 + 2 * a * rho; 
beta = sqrt(beta_2);
k_2 = 1 - alpha_2./beta_2;  
C = mu * I * pi^(-1);

% elliptic integrals of first and second kind for argument k^2
[K, E] = ellipke(k_2); 

Bx = C * x .* z .* (((a^2 + r.^2) .* E) - (alpha_2 .* K)) .* (2 * alpha_2 .* beta .* rho.^2).^(-1);

By = (y .* Bx) .* x.^(-1);

Bz = C * (((a^2 - r.^2) .* E) + (alpha_2 .* K)) .* (2 .* alpha_2 .* beta).^(-1);
end
