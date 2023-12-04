clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 14;
textFontSize = 10;
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
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);

tableFilename = 'ArticlePlot/stagnationDikeData.xlsx';
T = readtable(tableFilename);

sim_id = table2array(T(:, 1));
k  = table2array(T(:, 2));
ch20_sat = table2array(T(:, 6));
front = table2array(T(:, 4)) * 1e-3;
time = table2array(T(:, 7));
%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 89 / inch;
figY = 0.618 * figX;



fig = figure('Name', 'Contourf (front by k and ch2o)', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);

xx = k;
yy = ch20_sat;
zz = front;
F = scatteredInterpolant(xx, yy, zz);

[X, Y] = meshgrid(linspace(min(k), max(k), 7), linspace(min(ch20_sat), max(ch20_sat), 7));
Z = F(X, Y);
[C,h] = contourf(ax, X, Y, Z);
ylim(ax, [min(ch20_sat), max(ch20_sat)])
colormap(ax, 'parula')
h_cb = colorbar(ax);

xlabel(ax, "$k$, J/(m$^2$ K s)")
ylabel(ax, "$c_{H_2O}$, \%")
ylabel(h_cb, "$x$, km", FontSize=stFontSize, Interpreter='latex', Visible='on');
set(h_cb, TickLabelInterpreter='latex', FontSize=stFontSize);

hold on;
xx = k;
yy = time;
zz = ch20_sat;
F = scatteredInterpolant(xx, yy, zz);

T = [2 * 24 * 3600, 7 * 24 * 3600, 30 * 24 * 3600];
timeText = ["~2 days~", "~1 week~", "~1 month~"];
[X, tt] = meshgrid(linspace(min(k), max(k), 7), T);
Y = F(X, tt);

X = X';
Y = Y';
Xp = [X(5, 1), X(4, 2), X(2,3)];
Yp = [Y(5, 1), Y(4, 2), Y(2,3)];

plot(ax, X, Y, LineWidth=2, LineStyle='--', Color='w')
for i=1:length(T)
    xp = Xp(i);
    yp = Yp(i);
    strText = num2str(timeText(i));
    text(ax, xp, yp, strText, FontSize=textFontSize, HorizontalAlignment='center', EdgeColor='k', BackgroundColor='white', Clipping='on', Margin=0.5);
end


save_path = fullfile(save_dir, 'contour_k_ch2o_time.pdf');
exportgraphics(fig, save_path, ContentType='vector');