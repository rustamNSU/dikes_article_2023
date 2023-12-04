clear all;
close all;
addpath(genpath(pwd));

%% Define simulation directories and load data
timestep = 67;
sim_dir = fullfile(pwd, 'Simulations', 'simID_6000');
sim_filename = @(i) sprintf('timestep_%05d.mat', i);
load(fullfile(sim_dir, sim_filename(timestep)), 'FractureData');
load(fullfile(sim_dir, 'reservoir'), 'reservoir');

%% Fracture parameters
frac_elem = FractureData.channel;
X         = 1e-3 * FractureData.mesh.xc(frac_elem);  % elements center [km]
Width     = FractureData.width(frac_elem);           % elements width [m]
Pressure  = 1e-6 * FractureData.pressure(frac_elem); % fluid (magma) pressure [MPa]
Density   = FractureData.rho(frac_elem);             % fluid (magma) density [kg/m^3]
Viscosity = FractureData.mu(frac_elem);              % fluid (magma) viscosity [Pa*s]
SigmaH    = 1e-6 * reservoir.sigmaH(frac_elem);      % surrounded host-rock minimum pressure
T         = FractureData.temperature(frac_elem);
ElasticPressure = Pressure - SigmaH;
Channels  = FractureData.channel; 


%% Create figure space and plot data
fracture_fig = figure('Name', 'Fracture parameters', 'Units','normalized', 'Position',[.15,.15,.7,.9]);
sgtitle(fracture_fig, sprintf('Time = %ds,   time = %0.1fh', FractureData.time, FractureData.time / 3600));

%% Width
subplotWidth = subplot(3,2,1);
line_width = plot(subplotWidth, X, Width, 'linewidth', 2);
xlabel('Depth [km]');
ylabel('Width [m]');
grid on;
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title(subplotWidth, 'Width');

%% Pressure
subplotPressure = subplot(3,2,2);
line_pressure = plot(subplotPressure, X, Pressure, 'linewidth', 2);
xlabel('Depth [km]');
ylabel('Pressure [MPa]');
ylim([-1, Inf]);
grid on;
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title(subplotPressure, 'Pressure');

%% Elastic pressure
subplotElasticPressure = subplot(3,2,3);
line_elastic_pressure = plot(subplotElasticPressure, X, ElasticPressure, 'linewidth', 2);
xlabel('Depth [km]');
ylabel('Elastic pressure [MPa]');
ylim([-1, Inf]);
grid on;
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title('Elastic pressure');

%% Density
subplotRho = subplot(3,2,4);
line_rho = plot(subplotRho, X, Density, 'linewidth', 2);
xlabel('Depth [km]');
ylabel('Density [kg/m^3]');
grid on;
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title('Density');


%% Temperature
subplotT = subplot(3,2,5);
line_T = plot(subplotT, X, T, 'linewidth', 2);
xlabel('Depth [km]');
ylabel('T [C]');
grid on;
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title('Temperature');

%% Viscosity
subplotMu = subplot(3,2,6);
line_mu = semilogy(subplotMu, X, Viscosity, 'linewidth', 2);
grid on;
xlabel('Depth [km]');
ylabel('Viscosity [Pa\bullets]');
set(gca,'FontSize',12,'Fontweight', 'bold', 'linewidth',1.5);
title('Viscosity');