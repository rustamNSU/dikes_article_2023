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
sim_dir = {"simID_10121", "simID_10121"};
legendText = "";
saveFilename = 'gas_filtration';
timesteps = {0:1:370, 0:1:178};
sim_path = fullfile(sim_root_dir, sim_dir);

data = cellfun(@parameters_in_front, sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);


%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 183 / inch;
figY = 0.618 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);

ax1 = axes();
X1 = front(2:end) - front(end);
Y1 = Vg(2:end) ./ v + 1.0;
plot(ax1, X1, Y1, Color='black', LineStyle='-', Marker='o', MarkerIndices=1:1:length(X1));
ylabel(ax1, "$v_g / v$")
xlabel(ax1, '$x-x_s$, m')
set(ax1, xlim=[X1(end) - 15, X1(end)], ylim=[1, 10],...
    YScale='log', YTick=[1, 10, 100], YTickLabel=[1, 10, 100]);
grid on;



ax2 = axes('Position',[0.25 0.4 0.45 0.45]);
X2 = X1 * 1e-3;
S = 10;

yyaxis(ax2,'left');
    Y2l = kg(2:end);
    plot(ax2, X2, Y2l, Color='b', LineStyle='-');
    ylabel(ax2,"$k$, m$^2$");
    set(ax2, 'YColor','b', xlim=[X2(end) - S, X2(end)], YScale='log',...
        YTick=[1e-18, 1e-17, 1e-16, 1e-15, 1e-14]);
   
yyaxis(ax2,'right');
    Y2r = alpha(2:end);
    plot(ax2, X2, Y2r, Color='r', LineStyle='-');
    ylabel(ax2, "$\alpha$");
    set(ax2, 'YColor','r', xlim=[X2(end) - S, X2(end)]);

annotation('textarrow', [0.5, 0.7], [0.75, 0.71], 'String',{"second boiling at", "stagnation point"},...
    'LineWidth',1.5, 'HeadWidth',5, 'HeadLength',5,...
    'TextEdgeColor','black', 'interpreter','latex');
xlabel(ax2, '$x-x_s$, km')
set(ax2, 'FontSize', 14);

save_path = fullfile(save_dir, [sim_dir '_co2_sat.pdf']);
exportgraphics(fig, save_path);