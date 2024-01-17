function K = alldredge_inv_iter_first(A, ref_spectrum)
ref_spectrum = reshape(ref_spectrum,[],1);
K = (A'*A)\(A'*ref_spectrum);
end