function mat = gauss2d(max_val, mat, sigma, center)
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(max_val, R, C, sigma, center);