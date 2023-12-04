clear all;
close all;
clc;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 12;
stLegendFontSize = 10;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', stLegendFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

lineStyle = ["-", "-", "-", "-"];
lineMarkers = ["none", "none", "none", "none", "none", "none", "none", "none"];
xLim = [0.1, 100];
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images');
mkdir(save_dir);

simDirs = ["simID_1", "simID_2", "simID_3"];
timesteps = {1:5:3000, 1:5:1900, 1:5:3000};
legendText = [
    "csd",...
    "betaeq",...
    "linear, $tau=7$ days"];
saveFilename = 'compare_csd.pdf';
lineMarkers = {'none', 'none', 'none', 'none', 'none'};
simColors = {'k', 'b', 'g', 'r', 'y'};
xLim = [0.1, 80];
stLegendFontSize = 9;

sim_paths = fullfile(sim_root_dir, simDirs);
[front, time] = cellfun(@generate_front_depth_by_time, sim_paths, timesteps, 'UniformOutput', false);
[v, tv] = cellfun(@numeric_derivative, front, time, 'UniformOutput', false);

%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 150 / inch;
figY = 0.618 * figX;


fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);
yyaxis(ax,'left');
h_front = {};
for i = 1:length(simDirs)
    X = time{i} / 3600;
    Y = front{i} * 1e-3;
    h_1 = semilogx(ax, X, Y, Color=simColors{i}, LineStyle='-', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    h_front = [h_front, h_1];
    hold on;
end
ylim([-Inf, 0]);
set(ax, 'YColor','black');
ylabel(ax,"Depth, $l_t$ (km)", FontSize=stFontSize);

yyaxis(ax,'right');
for i = 1:length(simDirs)
    X = tv{i} / 3600;
%     X = 1:length(tv{i});
    Y = movmean(v{i},5);
    h_2 = loglog(ax, X, Y, Color=simColors{i}, LineStyle='--', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    hold on;
end
set(ax, 'xscale','log', 'yscale','log');
set(ax, 'YColor','black', 'YTick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],'YTickLabel',{'$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
ylim([1e-5, 10]);
xlim(xLim);
ylabel(ax,"Velocity (m/s)", FontSize=stFontSize);
ax.XGrid = 'on';
xlabel(ax, 'Time (h)', FontSize=stFontSize)
set(ax, 'XTick',[1, 10, 100],'XTickLabel',{'$1$', '$10$', '$100$'});

% hStop = xline(ax, 10000/3600, Color='r', LineStyle='--', LineWidth=2);
legend(h_front, legendText, Location='west', FontSize=stLegendFontSize);

% save_path = fullfile(save_dir, saveFilename);
% exportgraphics(fig, save_path);