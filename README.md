# kuslits-et-al-current-loops
Methodology searching for concentrations of electric currents in Earth's outer core using an alternate model of the geomagnetic field and deep learning.
Codes are developed and can run in a MATLAB R2014b environment.

# Alldredge-method
See the corresponding files in 'Alldredge-method/inversion'.\\
Compares our methodology against a more standard linear-non-linear inversion approach used by Alldredge (1987). \\
You can run the main script using 'inv_Alldredge.m'

# training-sample-generation
To generate a training set extract the file "test_core_new.zip" containing the historical GUFM-1 field model data and use "run_training_set_generation.m".

# GA-estimation
You can run a Genetic Algorithm based (GA) estimation for obtaining final estimations of loop parameters with "launch_and_prepare_ga_new.m"

# final-model-evaluation
Calculates the quality parameters of the final estimated loop model against the actual (ground truth) model using the script "est_vs_true.m".
