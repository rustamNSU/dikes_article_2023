clear all;
close all;
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
set(0, 'defaultLegendFontSize', stFontSize)
set(0, 'defaultLegendInterpreter', 'latex')


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_140001';
timestep = 450;
xLim = [-35, -5];

markers = {'o', '^'};
markers = {'none', 'none'};

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, generate_filename(timestep)), 'FractureData');
load(fullfile(sim_path, 'reservoir'), 'reservoir');



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 200 / inch;
figY = 0.618 * figX;

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
ElasticPressure = 1e-6 * reservoir.elasticityMatrix * Width';
alpha     = FractureData.alpha;
beta      = FractureData.beta;


%% Width and overpressure
fig = figure('Name', 'all', 'Units','inches', 'Position',[5, 2, figX, figY]);
sgtitle(fig,  sprintf('Time = %0.1fh', FractureData.time / 3600));
ax = subplot(2,2,1);
yyaxis(ax,'left');
    X = xc(ext_elem);
    Y = Width(ext_elem);
    h_w = plot(ax, X, Y, Color='b', LineStyle='-', Marker=markers{1}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Width (m)", FontSize=stFontSize);
    set(ax, 'YTick',[0, 2, 4, 6, 8, 10], 'YTickLabel',{'$1$', '$2$', '$4$', '$6$', '$8$', '$10$'});

yyaxis(ax,'right');
    X = xc(channel);
    Y = ElasticPressure(channel);
    h_ep = plot(ax, X, Y, Color='r', LineStyle='-', Marker=markers{2}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Overpressure (MPa)", FontSize=stFontSize);
    ylim([-10, 30])

ax.XGrid = 'on';
xlabel(ax, 'Depth (km)', FontSize=stFontSize)
legend([h_w, h_ep], "$w$", "$p_e$", Location='best')
set(ax, xlim=xLim)
title(ax, "\bf (a)", FontSize=stFontSize);


%% Density and buoyancy
ax = subplot(2,2,2);
yyaxis(ax,'left');
    X = xc(frac_elem);
    Y = Density(frac_elem);
    h_rho = plot(ax, X, Y, Color='b', LineStyle='-', Marker=markers{1}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Density (kg/m$^3$)", FontSize=stFontSize);

yyaxis(ax,'right');
    X = xc(frac_elem);
    Y = reservoirDensity(1e3*X);
    Y = Y - Density(frac_elem);
    h_p = plot(ax, X, Y, Color='r', LineStyle='-', Marker=markers{2}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Buoyancy (kg/m$^3$)", FontSize=stFontSize);

ax.XGrid = 'on';
xlabel(ax, 'Depth (km)', FontSize=stFontSize)
legend([h_rho, h_p], "$\rho$", "$\rho_{lit} - \rho$", Location='best')
set(ax, xlim=xLim)
title(ax, "\bf (b)", FontSize=stFontSize);

%% Temperature and viscosity
ax = subplot(2,2,3);
yyaxis(ax,'left');
    X = xc(frac_elem);
    Y = T(frac_elem);
    h_1 = plot(ax, X, Y, Color='b', LineStyle='-', Marker=markers{1}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Temperature (C$^\circ$)", FontSize=stFontSize);

yyaxis(ax,'right');
    X = xc(frac_elem);
    Y = Viscosity(frac_elem);
    h_2 = semilogy(ax, X, Y, Color='r', LineStyle='-', Marker=markers{2}, MarkerIndices=1:40:length(X));
    set(ax, 'YColor','black');
    ylabel(ax,"Viscosity (Pa$\cdot$s)", FontSize=stFontSize);
    set(ax, 'YTick',[1, 1e3, 1e6, 1e9, 1e12, 1e15], 'YTickLabel',{'$1$', '$10^{3}$', '$10^{6}$', '$10^{9}$', '$10^{12}$', '$10^{15}$'});
    ylim([1e3, 1e15])
%     set(ax, 'YTick',[1, 1e4, 1e8, 1e12, 1e16, 1e20],'YTickLabel',{'$1$', '$10^{4}$', '$10^{8}$', '$10^{12}$', '$10^{16}$', '$10^{20}$'});

ax.XGrid = 'on';
xlabel(ax, 'Depth (km)', FontSize=stFontSize)
legend([h_1, h_2], "$T$", "$\mu$", Location='west')
set(ax, xlim=xLim)
title(ax, "\bf (c)", FontSize=stFontSize);

%% alpha and beta
ax = subplot(2,2,4);
X = xc(frac_elem);
Y = alpha(frac_elem);
h_1 = plot(ax, X, Y, Color='b', LineStyle='-', Marker=markers{1}, MarkerIndices=1:40:length(X));
set(ax, 'YColor','black');
hold on;

X = xc(frac_elem);
Y = beta(frac_elem);
h_2 = plot(ax, X, Y, Color='r', LineStyle='-', Marker=markers{2}, MarkerIndices=1:40:length(X));

ax.XGrid = 'on';
ax.YGrid = 'on';
xlabel(ax, 'Depth (km)', FontSize=stFontSize)
legend([h_1, h_2], "$\alpha$", "$\beta$", Location='best')
ylabel(ax,"Volume concentrations", FontSize=stFontSize);
set(ax, xlim=xLim)
set(ax, 'YTick',[0, 0.2, 0.4, 0.6, 0.8, 1]);
title(ax, "\bf (d)", FontSize=stFontSize);


save_path = fullfile(save_dir, [sim_dir '_all.pdf']);
exportgraphics(fig, save_path);


% get density of reservoir (z is a depth in [m] with minus sign)
function rho = reservoirDensity(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    d = 1724.7;
    alpha = 0.001918;
    rho = rho1 - alpha * z - (rho1 - rho0) * exp(z / d);
end