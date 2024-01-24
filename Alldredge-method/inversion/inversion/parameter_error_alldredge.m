function MPE = parameter_error_alldredge(corr_pars, est_pars)
%%Calculates the relative parameter errors for a given estimation of the loop parameters 
%%using the Alldredge method
%% Finding the physically nearest true loop is not needed here as initial loops were fixed to be the same during the stability analysis
MPE = abs(corr_pars-est_pars)./abs(corr_pars);

end
