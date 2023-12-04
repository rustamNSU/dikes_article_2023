clear all;
close all;
clc;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
global stFontSize stLegendFontSize
stFontSize = 12;
stLegendFontSize = 10;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', stFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

sim_root_dir = fullfile(pwd, 'Simulations');
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 250 / inch;
figY = 0.618 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);

%% Conductivity
ax1 = subplot(2,2,1);
simDirs = ["simID_140011", "simID_140001", "simID_140010"];
timesteps = {1:3:4000, 1:2:1500, 1:2:1500};
legendText = [
    "$c_{H_2O}\!=\!2.35\%$", ...
    "$c_{H_2O}\!=\!6.24\%$", ...
    "$c_{H_2O}\!=\!8.66\%$"];
lineMarkers = {'none', 'none', 'none', 'none'};
simColors = {'b', 'k', 'r'};
xLim = [0.04, 200];
plotDikeFronts(ax1, simDirs, timesteps, legendText, ...
    lineMarkers, simColors, [-Inf, 0], [1e-5, 5], xLim);
title(ax1, "\bf (a)", FontSize=stFontSize);


%% ch2o
ax2 = subplot(2,2,2);
simDirs = ["simID_140012", "simID_140001", "simID_140013"];
timesteps = {1:2:1500, 1:2:1500, 1:2:1500};
legendText = [
    "$M = 5\times10^7$ kg$/$m", ...
    "$M = 10\times10^7$ kg$/$m", ...
    "$M = 15\times10^7$ kg$/$m"];
lineMarkers = {'none', 'none', 'none', 'none'};
simColors = {'b', 'k', 'r'};
xLim = [0.04, 100];
plotDikeFronts(ax2, simDirs, timesteps, legendText, ...
    lineMarkers, simColors, [-Inf, 0], [1e-5, 5], xLim);
title(ax2, "\bf (b)", FontSize=stFontSize);


%% Mass
ax3 = subplot(2,2,3);
simDirs = ["simID_140001", "simID_140014"];
timesteps = {1:2:1500, 1:2:1500};
legendText = ["$T_{ch}\!=\!850$ C$^\circ$", "$T_{ch}\!=\!900$ C$^\circ$"];
lineMarkers = {'none', 'none', 'none', 'none'};
simColors = {'k', 'r'};
xLim = [0.04, 100];
plotDikeFronts(ax3, simDirs, timesteps, legendText, ...
    lineMarkers, simColors, [-Inf, 0], [1e-5, 5], xLim);
title(ax3, "\bf (c)", FontSize=stFontSize);


%% Chamber temperature
ax4 = subplot(2,2,4);
simDirs = ["simID_140015", "simID_140001", "simID_140017"];
timesteps = {1:2:1500, 1:2:1500, 1:2:1500};
legendText = ["$t_{inj}\!=\!20000$ s", "$t_{inj}\!=\!10000$ s", "$t_{inj}\!=\!5000$ s"];
lineMarkers = {'none', 'none', 'none', 'none'};
simColors = {'b', 'k', 'r'};
xLim = [0.1, 100];
plotDikeFronts(ax4, simDirs, timesteps, legendText, ...
    lineMarkers, simColors, [-Inf, 0], [1e-5, 5], xLim);
title(ax4, "\bf (d)", FontSize=stFontSize);


save_path = fullfile(save_dir, "front_subplots_k_ch2o_M_Tch.pdf");
exportgraphics(fig, save_path);


function plotDikeFronts(ax, ...
    simDirs, ...
    timesteps, ...
    legendText, ...
    lineMarkers, ...
    simColors, ...
    ylimFront, ...
    ylimVelocity, ...
    xLim)
    global stFontSize stLegendFontSize
    sim_root_dir = fullfile(pwd, 'Simulations');
    sim_paths = fullfile(sim_root_dir, simDirs);
    [front, time] = cellfun(@generate_front_depth_by_time, sim_paths, timesteps, 'UniformOutput', false);
    [v, tv] = cellfun(@numeric_derivative, front, time, 'UniformOutput', false);

    yyaxis(ax,'left');
    h_front = {};
    for i = 1:length(simDirs)
        X = time{i} / 3600;
        Y = front{i} * 1e-3;
        h_1 = semilogx(ax, X, Y, Color=simColors{i}, LineStyle='-', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
        h_front = [h_front, h_1];
        hold on;
    end
    ylim(ylimFront);
    set(ax, 'YColor','black');
    ylabel(ax,"Depth, $l_t$ (km)", FontSize=stFontSize);

    yyaxis(ax,'right');
    for i = 1:length(simDirs)
        X = tv{i} / 3600;
        Y = movmean(v{i},5);
        h_2 = loglog(ax, X, Y, Color=simColors{i}, LineStyle='--', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
        hold on;
    end
    set(ax, 'xscale','log', 'yscale','log');
    set(ax, 'YColor','black', 'YTick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10],'YTickLabel',{'$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$1$', '$10$'});
    ylim(ylimVelocity);
    xlim(xLim);
    ylabel(ax,"Velocity (m/s)", FontSize=stFontSize);
    ax.XGrid = 'on';
    xlabel(ax, 'Time (h)', FontSize=stFontSize)
    set(ax, 'XTick',[0.1, 1, 10, 100],'XTickLabel',{'$0.1$', '$1$', '$10$', '$100$'});
    legend(h_front, legendText, Location='west', FontSize=stLegendFontSize);
end