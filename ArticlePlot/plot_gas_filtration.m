clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 22;
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
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_10121';
legendText = "";
saveFilename = 'gas_filtration';
timesteps = 0:1:95;
sim_path = fullfile(sim_root_dir, sim_dir);

[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[Vg, timeg, dP, rhog, alpha, kg] = gas_filtration_survey(sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);


%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 183 / inch;
figY = 0.618 * figX;

fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);
yyaxis(ax,'left');
    X = tv / 3600;
    Y = movmean(v,5);
    h_1 = loglog(ax, X, Y, Color='black', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    hold(ax, 'on');

    X = timeg(2:end) / 3600;
    Y = Vg(2:end);
    h_2 = loglog(ax, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    
    set(ax, 'YColor','black', 'YTick',[1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1],'YTickLabel',{'$10^{-15}$', '$10^{-12}$', '$10^{-9}$', '$10^{-6}$', '$10^{-3}$', '$10^{0}$'});
    ylabel(ax,"$v$, m/s", FontSize=stFontSize);
    set(ax, 'YColor','black');

yyaxis(ax,'right');
    X = timeg(2:end) / 3600;
    Y = kg(2:end) * 1.01324997e+15;
    h_3 = semilogx(ax, X, Y, Color='b', LineStyle='-', Marker='s', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"$k_m$, md", FontSize=stFontSize);

grid on;
xlabel(ax, '$t$, h', FontSize=stFontSize)
set(ax, 'XTick',[1, 10, 100],'XTickLabel',{'$1$', '$10$', '$100$'});
legend([h_1, h_2, h_3], "$v$", "$v_g$", "$k_m$", Location='west')

save_path = fullfile(save_dir, [sim_dir '_gas_filtration.pdf']);
exportgraphics(fig, save_path);