clear all;
close all;
addpath(genpath(pwd));

sim_root_dir = fullfile(pwd, 'Simulations');
simDirs = ["simID_16", "simID_17"];
timesteps = {0:3000, 0:3000};

sim_paths = fullfile(sim_root_dir, simDirs);
for i = 1:length(simDirs)
    sim_dir = sim_paths(i);
    timestep = timesteps{i};
    save_crystallization_data_by_time(sim_dir, timestep);
end