# kuslits-et-al-current-loops
Methodology searching for concentrations of electric currents in Earth's outer core using an alternate model of the geomagnetic field and deep learning.
Codes are developed and can run in MATLAB R2014b and python 3.10 environments.

## Alldredge-method
A more standard linear-non-linear inversion approach used by Alldredge (1987).
See the corresponding files in 'Alldredge-method/inversion'.
You can run the main script using 'inv_Alldredge.m'

## training-sample-generation
To generate a training set extract the file "test_core_new.zip" containing the historical GUFM-1 field model data and use 'run_training_set_generation.m'.

## MagneticDann (Deep Learning based estimation using a DANN trained UNet++ neural net)
To train the UNet++ neural network for obtaining estimated distributions loop parameters (as output maps) on the generated training set data use 'DannUnet.py'.

## loop-param-ranges
Calculation for deriving the range of potential loop parameters in the models. To perform this estimation run 'compute_par_ranges.m'.

## GA-estimation
You can run a Genetic Algorithm based (GA) estimation for obtaining final estimations of loop parameters with 'launch_and_prepare_ga_new.m'.

## final-model-evaluation
Calculates the quality parameters of the final estimated loop model against the actual (ground truth) model using the script 'est_vs_true.m'.

## tests
'eval_test_data_loop_discovery.m': Calculates the accuracy of the deep learning methodology for recovering initial loop positions on a given test set stored in 'current-loops-data.zip', and compares that against the recovery method of initial loop positions according to Alldredge (1987).
'eval_test_data_cross_correlations.m': Calculates the 2d Cross Correlation Coefficiens between true and recovered loop parameter distributions on a given test set stored in 'current-loops-data.zip'.
