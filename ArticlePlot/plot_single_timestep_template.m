clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', 14)
set(0, 'defaultLegendInterpreter', 'latex')


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


%% Base case
sim_dir = 'simID_130003';
legendText = "k=3 [J/(m^2*K)], X_{H_2O}=6.24%, T_0 = 850C\circ";
timestep = 1451;
xLim = [-35, 0];

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, generate_filename(timestep)), 'FractureData');
load(fullfile(sim_path, 'reservoir'), 'reservoir');



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 183 / inch;
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
T         = FractureData.temperature;
ElasticPressure = Pressure - SigmaH;


%% Width and overpressure
fig1 = figure('Name', 'width and overpressure', 'Units','inches', 'Position',[5, 2, figX, figY]);
ax1 = subplot(1,1,1);
yyaxis(ax1,'left');
    X = xc(ext_elem);
    Y = Width(ext_elem);
    h_w = plot(ax1, X, Y, Color='b', LineStyle='-', Marker='o', MarkerIndices=1:10:length(X));
    set(ax1, 'YColor','black');
    ylabel(ax1,"$w$, m", FontSize=14);

yyaxis(ax1,'right');
    X = xc(channel);
    Y = ElasticPressure(channel);
    h_ep = plot(ax1, X, Y, Color='r', LineStyle='-', Marker='^', MarkerIndices=1:10:length(X));
    set(ax1, 'YColor','black');
    ylabel(ax1,"$p_e$, MPa", FontSize=14);

grid on;
xlabel(ax1, '$x$, km', FontSize=14)
legend([h_w, h_ep], "$w$", "$p_e$", Location='Best')
set(ax1, xlim=xLim)

% save_path = fullfile(save_dir, [sim_dir '_w_overpressure.pdf']);
% exportgraphics(ax1, save_path);