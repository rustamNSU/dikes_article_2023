clear all;
close all;
clc;
addpath(genpath(pwd));

sim_root_dir = fullfile(pwd, 'Simulations');
kIndex = 0:5;
satIndex = [1, 2, 3, 4, 6];
satValue = [2.35, 4.12, 5.01, 6.24, 8.66];

[kGrid, satGrid] = ndgrid(kIndex, satIndex);
[~, satValueGrid] = ndgrid(kIndex, satValue);
simDirs = "simID_3" + string(kGrid) + string(satGrid) + "0";
simPaths = fullfile(sim_root_dir, simDirs);
totalFiles = arrayfun(@(folder) dir(folder), (fullfile(simPaths, '*.mat')), UniformOutput=false);
totalTimesteps = arrayfun(@(x) length(x{1}) - 2, totalFiles);
timesteps = arrayfun(@(x) 0:1:x-1, totalTimesteps, UniformOutput=false);

% [Vg, timeg, dP, rhog, alpha, kg] = cellfun(@gas_filtration_survey, sim_dirs, timesteps, 'UniformOutput', false);
[w, p, rho, mu, T, alpha, beta, cco2, ch2o, xh2o, xh2od, e, xFront, time, mcGas, mcCrystal] = cellfun(@get_data_from_front, simDirs, timesteps, 'UniformOutput', false);

stagnationDepth = arrayfun(@(x) max(x{1}), xFront);
velocity = cellfun(@numeric_derivative, xFront, time, UniformOutput=false);
stagnationTimeIndex = arrayfun(@(v) sum(v{1} > 1e-4) + 1, velocity);
stagnationTime = arrayfun(@(t, i) t{1}(i), time, stagnationTimeIndex);

tableFilename = 'ArticlePlot/stagnationDikeData.xlsx';
sim0 = reshape(simDirs, [numel(kGrid), 1]);
k1 = reshape(kGrid, [numel(kGrid), 1]);
sat2 = reshape(satGrid, [numel(kGrid), 1]);
front3 = reshape(stagnationDepth, [numel(kGrid), 1]);
time4 = reshape(stagnationTimeIndex, [numel(kGrid), 1]);
satValue5 = reshape(satValueGrid, [numel(kGrid), 1]);
time6 = reshape(stagnationTime, [numel(kGrid), 1]);

tableData = table(sim0, k1, sat2, front3, time4, satValue5, time6,...
    'VariableNames',["sim_ID", "k", "ch2o", "front", "stagnation time index (v<1e-3)", "water saturation, %", "stagnation time index (v<1e-3), s"]);
writetable(tableData,tableFilename, 'FileType','spreadsheet');


