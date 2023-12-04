clear all;
close all;
clc;
addpath(genpath(pwd));

sim_root_dir = fullfile(pwd, 'Simulations');
sim_dirs = {'simID_40100', 'simID_40204', 'simID_40203', 'simID_40202', 'simID_40111', 'simID_40200', 'simID_40201', 'simID_40112'};
page_names = ["equilibrium_crystallization", "tau=10000s", "tau=25000s", "tau=50000s", "tau=100000s", "tau=200000s", "tau=400000s", "tau=604800s"];
timesteps = {0:1:610, 0:1:610, 0:1:610, 0:1:610, 0:1:610, 0:1:610, 0:1:610, 0:1:610};
sim_paths = fullfile(sim_root_dir, sim_dirs);

front_data = cellfun(@parameters_in_front, sim_paths, timesteps, 'UniformOutput', false);


tableFilename = 'nonequilibrium_crystallization_front.xlsx';
for i = 1:length(sim_dirs)
    data = front_data{i};
    tableData = table(data.time', data.alpha', data.beta', data.betaeq', data.temperature', data.pressure', data.rho', data.mu', data.ch2od', data.cco2d', data.xh2od', data.xh2o', data.Mc', data.Mm', data.Mg', data.front',...
    'VariableNames',{'Time[s]', 'alpha[]', 'beta[]', 'betaeq[]', 'Temperature[K]', 'Pressure[Pa]', 'TotalDensity[kg/m^3]', 'Viscosity[Pa*s]', 'ch2o_liq[]', 'cco2_liq[]', 'xh2o_liq[]', 'xh2o[]', 'crystal_mass_content[]', 'melt_mass_content[]', 'gas_mass_content[]','Front[m]'});
    writetable(tableData,tableFilename, 'FileType','spreadsheet', 'Sheet',page_names(i));
end

