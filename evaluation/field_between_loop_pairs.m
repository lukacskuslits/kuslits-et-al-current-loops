R = 1;
I = 1;
mu = 1;

d1 = 0:1:10;
% d2 = 10:1:0;

B1 = mu*I*R^2./(2*(d1.^2+R^2).^(3/2));
B2 = mu*I*R^2./(2*((10-d1).^2+R^2).^(3/2));

B1 + B2