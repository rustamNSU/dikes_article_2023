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
sim_dir = 'simID_4300';
legendText = "k=3 [J/(m^2*K)], X_{H_2O}=6.24%, T_0 = 850C\circ";
timestep = 71;
timesteps = 0:1:370;
xLim = [-17, 0];

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, generate_filename(timestep)), 'FractureData');
load(fullfile(sim_path, 'reservoir'), 'reservoir');



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 3 * 183 / inch;
figY = 0.618 * figX / 3 * 2;


% Fracture parameters
frac_elem = FractureData.get_fracture_elements();
channel = FractureData.channel;
ext_elem = [frac_elem(1)-1, frac_elem, frac_elem(end)+1];

xc        = 1e-3 * FractureData.mesh.xc;
Width     = FractureData.width;           % elements width [m]
Pressure  = 1e-6 * FractureData.pressure; % fluid (magma) pressure [MPa]
Density   = FractureData.rho;             % fluid (magma) density [kg/m^3]
Viscosity = FractureData.mu;              % fluid (magma) viscosity [Pa*s]
SigmaH    = 1e-6 * reservoir.sigmaH;      % surrounded host-rock minimum pressure
T         = FractureData.temperature - 273.15;
ElasticPressure = Pressure - SigmaH;
alpha     = FractureData.alphaFrac';
beta      = FractureData.betaFrac';

[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);


%% Width and overpressure
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
sgtitle(fig, sprintf('Time = %ds,   time = %0.1fh', FractureData.time, FractureData.time / 3600));
ax1 = subplot(2,3,1);
yyaxis(ax1,'left');
    X = xc(ext_elem);
    Y = Width(ext_elem);
    h_w = plot(ax1, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    set(ax1, 'YColor','black');
    ylabel(ax1,"$w$, m", FontSize=stFontSize);

yyaxis(ax1,'right');
    X = xc(channel);
    Y = ElasticPressure(channel);
    h_ep = plot(ax1, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    set(ax1, 'YColor','black');
    ylabel(ax1,"$p_e$, MPa", FontSize=stFontSize);

grid on;
xlabel(ax1, '$x$, km', FontSize=stFontSize)
legend([h_w, h_ep], "$w$", "$p_e$", Location='north')
set(ax1, xlim=xLim)


%% Pressure and density
ax2 = subplot(2,3,2);
yyaxis(ax2,'left');
    X = xc(frac_elem);
    Y = Density(frac_elem);
    h_rho = plot(ax2, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    set(ax2, 'YColor','black');
    ylabel(ax2,"$\rho$, kg/m$^3$", FontSize=stFontSize);

yyaxis(ax2,'right');
    X = xc(frac_elem);
    Y = Pressure(frac_elem);
    h_p = plot(ax2, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    set(ax2, 'YColor','black');
    ylabel(ax2,"$p$, MPa", FontSize=stFontSize);

grid on;
xlabel(ax2, '$x$, km', FontSize=stFontSize)
legend([h_rho, h_p], "$\rho$", "$p$", Location='north')
set(ax2, xlim=xLim)


%% Temperature and viscosity
ax = subplot(2,3,4);
yyaxis(ax,'left');
    X = xc(frac_elem);
    Y = T(frac_elem);
    h_1 = plot(ax, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"$T$, C$^\circ$", FontSize=stFontSize);

yyaxis(ax,'right');
    X = xc(frac_elem);
    Y = Viscosity(frac_elem);
    h_2 = semilogy(ax, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"$\mu$, Pa s", FontSize=stFontSize);

grid on;
xlabel(ax, '$x$, km', FontSize=stFontSize)
legend([h_1, h_2], "$T$", "$\mu$", Location='north')
set(ax, xlim=xLim)


%% alpha and beta
ax = subplot(2,3,5);
X = xc(frac_elem);
Y = alpha;
h_1 = plot(ax, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
set(ax, 'YColor','black');
hold on;

X = xc(frac_elem);
Y = beta;
h_2 = plot(ax, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));

grid on;
xlabel(ax, '$x$, km', FontSize=stFontSize)
legend([h_1, h_2], "$\alpha$", "$\beta$", Location='north')
set(ax, xlim=xLim)

%% fracture front
ax = subplot(2,3,3);
yyaxis(ax,'left');
    X = time(2:end) / 3600;
    Y = front(2:end) * 1e-3;
    h_1 = semilogx(ax, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"$l_t$, km", FontSize=stFontSize);

yyaxis(ax,'right');
    X = tv / 3600;
    Y = movmean(v,5);
    h_2 = loglog(ax, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    set(ax, 'YColor','black', 'YTick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],'YTickLabel',{'$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
    ylabel(ax,"$v$, m/s", FontSize=stFontSize);

grid on;
xlabel(ax, '$t$, h', FontSize=stFontSize)
legend([h_1, h_2], "$l_t$", "$v$", Location='west')

save_path = fullfile(save_dir, [sim_dir '_all_front_t_' int2str(timestep) '.pdf']);
exportgraphics(fig, save_path);