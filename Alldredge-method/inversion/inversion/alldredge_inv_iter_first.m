function K = alldredge_inv_iter_first(A, ref_spectrum)
%% Computes the linear inversion of normalized magnetic moment K 
%% instability  and divergence from the correct parameters
%% remains even when strong regularization is applied here (see Hoerl and Kennard 1970)
% Input params:
%--------------
% A: systems matrix of Alldredge's equation
% ref_spectrum: correct SHCs
% Ouptut values:
%---------------
% K: normalized magnetic moments

lam=1e-3; %Regularization factor
I = ones(length(A(1,:)));
ref_spectrum = reshape(ref_spectrum,[],1);
K = (A'*A+lam*I)\(A'*ref_spectrum);
end