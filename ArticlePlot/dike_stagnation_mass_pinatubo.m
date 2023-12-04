clear all;
close all;
clc;
addpath(genpath(pwd));

sim_root_dir = fullfile(pwd, 'Simulations');
MIndex = 1:5;
satIndex = [1, 2, 3, 4, 5, 6, 7];
satValue = [2.37, 3.85, 5.13, 6.23, 8.13, 9.62, 11.15];

[MGrid, satGrid] = ndgrid(MIndex, satIndex);
[~, satValueGrid] = ndgrid(MIndex, satValue);
simDirs = "simID_1430" + string(MGrid) + string(satGrid);
simPaths = fullfile(sim_root_dir, simDirs);
totalFiles = arrayfun(@(folder) dir(folder), (fullfile(simPaths, '*.mat')), UniformOutput=false);
totalTimesteps = arrayfun(@(x) length(x{1}) - 4, totalFiles);
timesteps = arrayfun(@(x) 0:3:x-1, totalTimesteps, UniformOutput=false);

% [Vg, timeg, dP, rhog, alpha, kg] = cellfun(@gas_filtration_survey, sim_dirs, timesteps, 'UniformOutput', false);
[xFront, time] = cellfun(@get_data_from_front, simDirs, timesteps, 'UniformOutput', false);

stagnationDepth = arrayfun(@(x) max(x{1}), xFront);
velocity = cellfun(@numeric_derivative, xFront, time, UniformOutput=false);
stagnationTimeIndex = arrayfun(@(v) sum(v{1} > 1e-4) + 1, velocity);
stagnationTime = arrayfun(@(t, i) t{1}(i), time, stagnationTimeIndex);

tableFilename = 'ArticlePlot/stagnationDikeData_ch2o_M_pinatuba.xlsx';
sim0 = reshape(simDirs, [numel(MGrid), 1]);
k1 = reshape(MGrid, [numel(MGrid), 1]);
sat2 = reshape(satGrid, [numel(MGrid), 1]);
front3 = reshape(stagnationDepth, [numel(MGrid), 1]);
time4 = reshape(stagnationTimeIndex, [numel(MGrid), 1]);
satValue5 = reshape(satValueGrid, [numel(MGrid), 1]);
time6 = reshape(stagnationTime, [numel(MGrid), 1]);

tableData = table(sim0, k1, sat2, front3, time4, satValue5, time6,...
    'VariableNames',["sim_ID", "k", "ch2o", "front", "stagnation time index (v<1e-4)", "water saturation, %", "stagnation time (v<1e-4), s"]);
writetable(tableData,tableFilename, 'FileType','spreadsheet');


function [xFront, time] = get_data_from_front(simDir, timesteps)
    %% Load files from folder
    sim_filename = @(i) sprintf('timestep_%05d.mat', i);
    xFront= zeros(1, length(timesteps));
    time  = zeros(1, length(timesteps));
    
    
    iter = 1;
    for i = timesteps
        filepath = fullfile(simDir, sim_filename(i));
        load(filepath, 'FractureData');
        xFront(iter)= FractureData.front(end);
        time(iter)  = FractureData.time;
        iter = iter + 1;
    end
end
