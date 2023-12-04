%% Run simulation of magma dikes propogation from input data
clear all;
close all;
clc;
addpath(genpath(pwd));
warning off;

%% Load input data (pre-saved data)
% inputPath = {
%     'input/input_simID_143018.mat', 'input/input_simID_143028.mat', 'input/input_simID_143038.mat', 'input/input_simID_143048.mat', 'input/input_simID_143058.mat', ...
%     'input/input_simID_143019.mat', 'input/input_simID_143029.mat', 'input/input_simID_143039.mat', 'input/input_simID_143049.mat', 'input/input_simID_143059.mat'};
% saveFilename = {
%     'simID_143018', 'simID_143028', 'simID_143038', 'simID_143048', 'simID_143058',...
%     'simID_143019', 'simID_143029', 'simID_143039', 'simID_143049', 'simID_143059'};
inputPath = {'input/input_simID_10.mat', 'input/input_simID_18.mat'};
saveFilename = {'simID_10', 'simID_18'};
for i=1:length(inputPath)
    try
    	simulation_runner(inputPath{i}, saveFilename{i});
    catch 
    	
    end
end
 