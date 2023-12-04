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
save_dir = fullfile(pwd, 'images/Article2022-review');
mkdir(save_dir);

tableFilename = 'ArticlePlot/stagnationDikeData_ch2o_M_pinatuba.xlsx';
T = readtable(tableFilename);

mass = [2, 3, 4, 5, 6] * 2.5;
sim_id = table2array(T(:, 1));
M  = table2array(T(:, 2));
M = mass(M)';
ch20_sat = table2array(T(:, 6));
front = table2array(T(:, 4)) * 1e-3;
time = table2array(T(:, 7));
%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 200 / inch;
figY = 0.6 * figX;



fig = figure('Name', 'Contourf (front by M and ch2o)', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);

xx = M;
yy = ch20_sat;
zz = front;

xx = reshape(xx, 5, []);
yy = reshape(yy, 5, []);
zz = reshape(zz, 5, []);
F = griddedInterpolant(xx, yy, zz, 'cubic');

[X, Y] = meshgrid(linspace(min(M), max(M), 50), linspace(min(ch20_sat), max(ch20_sat), 50));
Z = F(X, Y);
[C,h] = contourf(ax, X, Y, Z, 50);
h.LineWidth = 0.01;
h.LineStyle = '-';
h.LineColor = 'none';
% h.ShowText = 'on';
% ylim(ax, [min(ch20_sat), max(ch20_sat)])
ylim(ax, [min(ch20_sat), 10]);
colormap(ax, 'turbo');
h_cb = colorbar(ax);

xlabel(ax, "Total magma mass $\times 10^{-7}$ (kg$/$m)")
ylabel(ax, "H$_2$O content (wt\%)")
ylabel(h_cb, "Depth (km)", FontSize=stFontSize, Interpreter='latex', Visible='on');
set(h_cb, TickLabelInterpreter='latex', FontSize=stFontSize);

hold on;
xx = M;
yy = time;
zz = ch20_sat;
F = scatteredInterpolant(xx, yy, zz);

T = [0.25 * 24 * 3600, 0.5 * 24 * 3600, 1 * 24 * 3600, 7 * 24 * 3600];
timeText = ["~$0.25$ day~", "~$0.5$ day~", "~1 day~", "~7 days~"];
[X, tt] = meshgrid(linspace(min(M), max(M), 7), T);
Y = F(X, tt);

X = X';
Y = Y';
% Xp = [X(5, 1)-1, (X(5, 1)+X(4, 2))/2-0.8, (X(4, 2)+X(2,3))/2-2.1];
% Yp = [Y(4, 1), Y(4, 2), Y(4, 3)+0.7];
Xp = [X(5,1), X(5, 2)-1, (X(5, 2)+X(4, 3))/2-0.8, (X(4, 3)+X(2,4))/2-2.1];
Yp = [Y(4,1)-0.3, Y(4, 2), Y(4, 3), Y(4, 4)+0.7];

plot(ax, X, Y, LineWidth=1.5, LineStyle='--', Color='k')
for i=1:length(T)
    xp = Xp(i);
    yp = Yp(i);
    strText = num2str(timeText(i));
    text(ax, xp, yp, strText, FontSize=textFontSize, HorizontalAlignment='center', EdgeColor='M', BackgroundColor='white', Clipping='on', Margin=0.5);
end


save_path = fullfile(save_dir, 'contour_M_ch2o_time.pdf');
exportgraphics(fig, save_path, 'Resolution',500);
% exportgraphics(fig, save_path, ContentType='vector');