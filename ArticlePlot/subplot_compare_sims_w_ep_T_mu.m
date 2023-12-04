clear all;
close all;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 22;
stFontSizeLog = 18;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', stFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

lineStyle = ["-", "--", ":", "-."];
lineMarkers = {'o', '^', 's', '*'};
simColors = {'b', 'black', 'r'};


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);


% %% Different k1c
% simDirs = ["simID_10121", "simID_31121"];
% legendText = ["$K_{1c}\!=\!1$ MPa$\cdot$m$^{1/2}$", "$K_{1c}\!=\!100$ MPa$\cdot$m$^{1/2}$"];
% saveFilename = 'compare_reservoir_toughness_parameters.pdf';
% simColors = {'black', 'r'};
% lineMarkers = {'^', 's'};
% timestep = 90;
% xLim = [-35, -5];

% %% Different E
% simDirs = ["simID_10121", "simID_32121"];
% legendText = ["$E\!=\!15$ GPa", "$E\!=\!45$ GPa"];
% saveFilename = 'compare_reservoir_young_parameters.pdf';
% simColors = {'black', 'r'};
% lineMarkers = {'^', 's'};
% timestep = 90;
% xLim = [-35, -5];





sim_paths = fullfile(sim_root_dir, simDirs);
sim_filenames = repmat({generate_filename(timestep)}, 1, length(sim_paths));
for i = 1:length(sim_filenames)
    sim_filenames{i} = fullfile(sim_paths{i}, sim_filenames{i});
end

%% Load compared files
FractureData = cell(1, length(sim_filenames));
for i = 1:length(sim_filenames)
    FractureData{i} = load(sim_filenames{i});
    FractureData{i} = FractureData{i}.FractureData;
end

reservoir = cellfun(@(simID) load(fullfile(sim_root_dir, simID, 'reservoir')), simDirs, 'UniformOutput', false);


%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 2 * 183 / inch;
figY = 0.618 * figX;

fig = figure('Name', 'all', 'Units','inches', 'Position',[5, 2, figX, figY]);
hLine1 = [];
hLine2 = [];
hLine3 = [];
hLine4 = [];
for i = 1:length(sim_filenames)
    % Fracture parameters
    frac_elem = FractureData{i}.get_fracture_elements();
    channel = FractureData{i}.channel;
    ext_elem = [frac_elem(1)-1, frac_elem, frac_elem(end)+1];

    xc        = 1e-3 * FractureData{i}.mesh.xc;
    Width     = FractureData{i}.width;           % elements width [m]
    Pressure  = 1e-6 * FractureData{i}.pressure; % fluid (magma) pressure [MPa]
    Density   = FractureData{i}.rho;             % fluid (magma) density [kg/m^3]
    Viscosity = FractureData{i}.mu;              % fluid (magma) viscosity [Pa*s]
    SigmaH    = 1e-6 * reservoir{i}.reservoir.sigmaH;      % surrounded host-rock minimum pressure
    T         = FractureData{i}.temperature - 273.15;
    ElasticPressure = Pressure - SigmaH;
    alpha     = FractureData{i}.alphaFrac';
    beta      = FractureData{i}.betaFrac';


    %% Width
    ax1 = subplot(2,2,1);
    X = xc(ext_elem);
    Y = Width(ext_elem);
    h_w = plot(ax1, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    xlabel(ax1, '$x$, km', FontSize=stFontSize);
    ylabel(ax1,"$w$, m", FontSize=stFontSize);
    grid on;
    set(ax1, xlim=xLim);
    title(ax1, "\bf (a)", FontSize=stFontSize);
    hold(ax1, 'on');


    %% Overpressure
    ax2 = subplot(2,2,2);
    X = xc(channel);
    Y = ElasticPressure(channel);
    h_ep = plot(ax2, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    grid on;
    xlabel(ax2, '$x$, km', FontSize=stFontSize);
    ylabel(ax2,"$p_e$, MPa", FontSize=stFontSize);
    set(ax2, xlim=xLim);
    title(ax2, "\bf (b)", FontSize=stFontSize);
    hold(ax2, 'on');

    %% Temperature
    ax3 = subplot(2,2,3);
    X = xc(frac_elem);
    Y = T(frac_elem);
    h_T = plot(ax3, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    grid on;
    xlabel(ax3, '$x$, km', FontSize=stFontSize);
    ylabel(ax3,"$T$, C$^\circ$", FontSize=stFontSize);
    set(ax3, xlim=xLim);
    title(ax3, "\bf (c)", FontSize=stFontSize);
    hold(ax3, 'on');

    %% Temperature
    ax4 = subplot(2,2,4);
    X = xc(frac_elem);
    Y = Viscosity(frac_elem);
    h_mu = semilogy(ax4, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    set(ax4, 'YTick',[1e4, 1e7, 1e10, 1e13, 1e16],'YTickLabel',{'$10^{4}$', '$10^{7}$', '$10^{10}$', '$10^{13}$', '$10^{16}$'});
    grid on;
    xlabel(ax4, '$x$, km', FontSize=stFontSize);
    ylabel(ax4,"$\mu$, Pa$\cdot$s", FontSize=stFontSize);
    set(ax4, xlim=xLim);
    title(ax4, "\bf (d)", FontSize=stFontSize);
    hold(ax4, 'on');
    
    hLine1 = [hLine1, h_w];
    hLine2 = [hLine2, h_ep];
    hLine3 = [hLine3, h_T];
    hLine4 = [hLine4, h_mu];
end
legend(ax1, hLine1, legendText, Location='best')
legend(ax2, hLine2, legendText, Location='best')
legend(ax3, hLine3, legendText, Location='best')
legend(ax4, hLine4, legendText, Location='best')
save_path = fullfile(save_dir, saveFilename);
exportgraphics(fig, save_path);


function rho = rhoLith(x_km)
    alpha=0.001918;
    d=1.7247;
    rho1=2.7323;
    rho0=2.0929;
    rho = 1e3 * (rho1 - alpha * x_km - (rho1 - rho0) * exp(x_km / d));
end