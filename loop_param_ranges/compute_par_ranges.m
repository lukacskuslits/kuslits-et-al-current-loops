% Regularizing the model with parameter bounds
clear all

% X = [1:Ilmin 2:Ilmax 3:Rlmin 4:Rlmax 5:nmin 6:nmax] 7:dmin 8:dmax]
%[]

GMD = 7.8e22;%8.1e22;
mu = 4*pi()*1e-7;
Roc = 3480 * 1e3;
Rdepth = 2680 * 1e3;%2680*1e3
Bmin = 1e-3;
Bmax = 1e-2;

norm1 = [1e24; 1e6; 1e6; 1e3; 1e3];
%norm1 = [1;1;1;1;1];

GMD = GMD/norm1(1);
Roc = Roc/norm1(2);
Rdepth = Rdepth/norm1(3);
Bmin = Bmin/norm1(4);
Bmax = Bmax/norm1(5);

Cw = (4/3)*pi()*(Roc^3-Rdepth^3)*(sqrt(2))^(-1);
C = Cw*sqrt(2);


bound_fn = @(X) [mu*X(1)*(2*X(4))^(-1)-Bmin; mu*X(2)*(2*X(3))^(-1)-Bmax; pi()*X(2)*X(4)^2*X(5)-GMD; pi()*X(1)*X(3)^2*X(6)-GMD;...
                 mu*X(4)^2*X(2)*(X(4)^2+X(5)^(-2/3)*Cw^(2/3))^(-3/2)-Bmin; mu*X(1)*(2*X(3))^(-1)+mu*X(1)*X(3)^2*(X(3)^2+X(6)^(-2/3)*C^(3/2))^(-3/2)-Bmax];
             
jacobi_fn = @(X) [mu*(2*X(4))^(-1), 0, 0, -mu*X(1)*(2*X(4)^2)^(-1), 0, 0;...
                  0, mu*(2*X(3))^(-1), -mu*X(2)*(2*X(3))^(-1), 0, 0, 0;...
                  0, pi()*X(4)^2*X(5), 0, 2*pi()*X(4)*X(2)*X(5), pi()*X(2)*X(4)^2, 0;...
                  pi()*X(3)^2*X(6), 0, 2*pi()*X(3)*X(1)*X(6), 0, 0, pi()*X(1)*X(3)^2;...
                  0, mu*X(4)^2*(X(4)^2+X(5)^(-2/3)*Cw^(2/3))^(-3/2), 0,...
                  mu*2*X(4)*X(2)*(X(4)^2+X(5)^(-2/3)*Cw^(2/3))^(-3/2)+mu*X(4)^2*X(2)*(-3/2)*(X(4)^2+X(5)^(-2/3)*Cw^(2/3))^(-5/2)*2*X(4),...
                  mu*X(4)^2*X(2)*(X(4)^2+X(5)^(-2/3)*Cw^(2/3))^(-5/2)*Cw^(2/3)*X(5)^(-5/3), 0;...
                  mu*(2*X(3))^(-1)+mu*X(3)^2*(X(3)^2+X(6)^(-2/3)*C^(3/2))^(-3/2), 0,...
                  (-1/2)*mu*X(1)*X(3)^(-2)+mu*2*X(3)*X(1)*(X(3)^2+X(6)^(-2/3)*C^(2/3))^(-3/2)+mu*X(3)^2*X(1)*(-3/2)*(X(3)^2+X(6)^(-2/3)*C^(2/3))^(-5/2)*2*X(3),...
                  0, 0, mu*X(3)^2*X(1)*(X(3)^2+X(6)^(-2/3)*C^(2/3))^(-5/2)*C^(2/3)*X(6)^(-5/3)];


X = [1e7; 1e9; 1e2; 1e6; 1e1; 1e3];
norm2 = [1e9; 1e9; 1e6; 1e6; 1e3; 1e3];
%norm2 = [1;1;1;1;1;1];
X = X./norm2;

no_itr = 1.5e3;
error = 2e-5;%3.5e-4; %3.344e-4;%
lam = 1;

%[point,no_itr,error_out]=NewtonRaphson_nl(X,bound_fn,jacobi_fn,no_itr,error,lam);
[point,no_itr,error_out]=NewtonRaphson_nl_print(X,bound_fn,jacobi_fn,no_itr,error,lam);
disp(point(1)*norm2(1))
disp(point(2)*norm2(2))
disp(point(3)*norm2(3))
disp(point(4)*norm2(4))
disp(point(5)*norm2(5))
disp(point(6)*norm2(6))
