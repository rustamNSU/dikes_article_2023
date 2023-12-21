
clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', 14)
set(0, 'defaultLegendInterpreter', 'latex')


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_5';
timesteps = [0:20:1000, 1040:40:3200];
xLim = [-400, 0];

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, 'reservoir'), 'reservoir');

[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);

sim_dir = 'simID_4';
timesteps = [0:20:1000, 1040:40:3200];
sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, 'reservoir'), 'reservoir');
[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[vp, tvp] = numeric_derivative(front, time);



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 183 / inch;
figY = 0.618 * figX;
% figX = 250 / inch;
% figY = 0.3*figX;


%% Width and overpressure
fig = figure('Name', 'width and overpressure', 'Units','inches', 'Position',[5, 2, figX, figY]);
ax2 = subplot(1,1,1);
X = tv;
Y = movmean(v,5);
h1 = loglog(ax2, X, Y, Color='b', LineStyle='-');
X = tvp;
Y = movmean(vp,5);
hold on;
h2 = loglog(ax2, X, Y, Color='r', LineStyle='-');
set(ax2, 'YColor','black');
grid on;
xlabel(ax2, 'Time (s)', FontSize=14);
ylabel(ax2,"Velocity (m/s)", FontSize=14);
set(ax2, xlim=[1e3, 2e5]);
fieldData = table2array(readtable('davis_num.csv'));
xx = fieldData(:,1);
yy = fieldData(:,2);
h3 = loglog(ax2, xx, yy, Color='green', LineStyle='--');

fieldData = table2array(readtable('davis_analyt.csv'));
xx = fieldData(:,1);
yy = fieldData(:,2);
h4 = loglog(ax2, xx, yy, Color='k', LineStyle='-.');
legend([h1, h2, h3, h4], ["this study, KGD", "this study, PKN", "num. (Davis, et. al. 2022)", "analyt. (Davis, et. al. 2022)"], Location='Best');
save_path = fullfile(save_dir, 'verification_Davis2022.pdf');
exportgraphics(fig, save_path);