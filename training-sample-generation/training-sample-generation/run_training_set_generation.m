
set = 1:5:6; %number of loops in the training set by batches - originally: 25:10:175; 
batch_size = 2; % size of each batch
attenuation = 1100; % factor (\gamma) for reducing the maximum potential rate of change in the loop currents (dI/dt)
deg_res = 2; % angular resolution of images (maps) in degrees
[sim_data, sample_numbers] = test_train_generator_new(set, batch_size, attenuation, deg_res);
