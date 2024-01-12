function MPE = parameter_error_alldredge(corr_pars, est_pars)

MPE = mean(mean(abs(corr_pars-est_pars)./abs(corr_pars)));

end