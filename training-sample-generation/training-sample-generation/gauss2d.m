function mat = gauss2d(max_val, mat, sigma, center)
%% Projects 2d Gaussian distribution values on a grid
%Iput params:
%------------
% max_val: value at the center (1 when data are normalized)
% mat: input matrix containing original {0,1} values for the grid
% sigma: standard deviation
% center: coordinate values of the center of the distribution (xc, yc)
%Output values:
%--------------
% mat: output matrix containing Gaussian distribution values

gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(max_val, R, C, sigma, center);
