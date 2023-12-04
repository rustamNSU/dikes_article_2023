clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 16;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', stFontSize)
set(0, 'defaultLegendInterpreter', 'latex')


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022-review');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_143001';
timesteps = 1:4:2000;
sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, 'reservoir'), 'reservoir');



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 150 / inch;
figY = 0.618 * figX;

[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);

fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);
yyaxis(ax,'left');
    X = time(2:end) / 3600;
    Y = front(2:end) * 1e-3;
    h_1 = semilogx(ax, X, Y, Color='b', LineStyle='-', Marker='none', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Depth, $l_t$ (km)", FontSize=stFontSize);
    ylim([-33, 0])

yyaxis(ax,'right');
    X = tv / 3600;
    Y = movmean(v,5);
    h_2 = loglog(ax, X, Y, Color='r', LineStyle='-', Marker='none', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black', 'YTick',[1e-6, 1e-4, 1e-2, 1, 10],'YTickLabel',{'$10^{-6}$', '$10^{-4}$', '$10^{-2}$', '1', '10'});
    ylabel(ax,"Velocity, $v$ (m/s)", FontSize=stFontSize);
    ylim([1e-5, 10]);

ax.XGrid = 'on';
xlim([0.1, 100]);
xlabel(ax, 'Time (h)', FontSize=stFontSize)
set(ax, 'XTick',[1, 10, 100],'XTickLabel',{'$1$', '$10$', '$100$'});
legend([h_1, h_2], "front depth", "front velocity", Location='west')

save_path = fullfile(save_dir, [sim_dir '_front_and_velocity.pdf']);
exportgraphics(fig, save_path);