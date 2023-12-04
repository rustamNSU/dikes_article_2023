clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 16;
textFontSize = 14;
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
save_dir = fullfile(pwd, 'images/bristol');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_130003';
legendText = sim_dir;
timestep = 1500;

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, generate_filename(timestep)), 'FractureData');
load(fullfile(sim_path, 'reservoir'), 'reservoir');



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figY = 183 / inch;
figX = 1 * figY;


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
Tmagma         = FractureData.temperature - 273.15;
ElasticPressure = Pressure - SigmaH;
Thost = FractureData.rockTemperature.T - 273.15;
khost = FractureData.rockTemperature.k;


%% Width and overpressure
fig = figure('Name', 'width and overpressure', 'Units','inches', 'Position',[5, 2, figX, figY]);
sgtitle(fig,  sprintf('Time = %0.1fh', FractureData.time / 3600), FontSize=stFontSize);
ax = subplot(1,5,[1,2]);
xx = xc;
yy = FractureData.rockTemperature.yc;
[X, Y] = meshgrid(xx, yy);
Z = Thost';
[C,h] = contourf(ax, Y, X, Z, LevelList=linspace(min(Z, [], 'all'), max(Z, [], 'all'), 50), EdgeColor='none');
colormap(ax, 'turbo')
h_cb = colorbar(ax);
xlabel(ax, "Wall rocks (horizontal) [m]")
ylabel(ax, "Depth [km]")
ylabel(h_cb, "Temperature, C$^\circ$", FontSize=stFontSize, Interpreter='latex', Visible='on');
set(h_cb, TickLabelInterpreter='latex', FontSize=stFontSize);
xlim([0, 1])

ax = subplot(1,5,[4,5]);
xx = xc;
yy = FractureData.rockTemperature.yc;
[X, Y] = meshgrid(xx, yy);
Z = khost';
[C,h] = contourf(ax, Y, X, Z, LevelList=linspace(min(Z, [], 'all'), max(Z, [], 'all'), 50), EdgeColor='none');
colormap(ax, 'turbo')
h_cb = colorbar(ax);
xlabel(ax, "Wall rocks (horizontal) [m]")
ylabel(ax, "Depth [km]")
ylabel(h_cb, "Thermal conductivity, W$/$(m$\cdot$K)", FontSize=stFontSize, Interpreter='latex', Visible='on');
set(h_cb, TickLabelInterpreter='latex', FontSize=stFontSize);
xlim([0, 1])


save_path = fullfile(save_dir, [sim_dir '_rock_temperature.pdf']);
exportgraphics(fig, save_path);
save_path = fullfile(save_dir, [sim_dir '_rock_temperature.jpg']);
exportgraphics(fig, save_path, 'Resolution',600);