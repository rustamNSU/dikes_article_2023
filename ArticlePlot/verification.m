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
sim_dir = 'simID_1';
legendText = "K=1";
timesteps_to_plot = [200, 600, 1000, 1660];
timesteps = 0:10:1660;
xLim = [-30, -10];

sim_path = fullfile(sim_root_dir, sim_dir);
load(fullfile(sim_path, 'reservoir'), 'reservoir');

[front, time] = generate_front_depth_by_time(sim_path, timesteps);
[v, tv] = numeric_derivative(front, time);



%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 183 / inch;
figY = 0.618 * figX;
figX = 250 / inch;
figY = 0.3*figX;

% Roper, Lister (2007)
fieldData = table2array(readtable('Roper_Lister_K=1_fig.2a.csv','DecimalSeparator',','));
% fieldData = table2array(readtable('2d_lister.csv','DecimalSeparator',','));
xx = fieldData(:,1);
yy = fieldData(:,2);


%% Width and overpressure
fig = figure('Name', 'width and overpressure', 'Units','inches', 'Position',[5, 2, figX, figY]);
ax1 = subplot(1,2,1);
xScale = 1825.74; % see verification_parameters.nb
X = -xx * xScale / 1e3 - 12.27;
Y = 2*yy;
X = [X; -40];
Y = [Y; 2];
h2 = plot(ax1, X, Y, Color='r', LineStyle='--');
hold on;
grid on;

% Fracture parameters
for i = timesteps_to_plot
    load(fullfile(sim_path, generate_filename(i)), 'FractureData');
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

    X = xc(ext_elem);
    Y = Width(ext_elem);
    h1 = plot(ax1, X, Y, Color='b', LineStyle='-');
    xp = X(end);
    yp = 1;
    time = FractureData.time;
    strText = sprintf('%0.2fh', FractureData.time / 3600);
    text(ax1, xp, yp, strText, FontSize=12, HorizontalAlignment='center', EdgeColor='k', BackgroundColor='white', Clipping='on', Margin=0.5);
    hold on;
end
set(ax1, 'YColor','black');
ylabel(ax1,"Width (m)", FontSize=14);

xlabel(ax1, 'Depth (km)', FontSize=14)
legend([h1, h2], "this study", "Roper, Lister (2007)", Location='Best')
set(ax1, xlim=xLim)
set(ax1, ylim=[0, 4])

ax2 = subplot(1,2,2);
X = tv / 3600;
Y = movmean(v,5);
h3 = plot(ax2, X, Y, Color='b', LineStyle='-');
set(ax2, 'YColor','black');
hold on;

c = 10;
h4 = yline(ax2, c, Color='r', LineStyle='--',  LineWidth=2);
legend([h3, h4], "this study", "Roper, Lister (2007)", Location='Best');
grid on;
xlabel(ax2, 'Time (h)', FontSize=14);
ylabel(ax2,"Velocity (m/s)", FontSize=14);
set(ax2, ylim=[7, 20]);

% save_path = fullfile(save_dir, 'verification_RoperLister2007.pdf');
% exportgraphics(fig, save_path);