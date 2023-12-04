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

lineStyle = ["-", "--", ":", "-."];
lineMarkers = {'o', '^', 's', '*'};
simColors = {'black', 'b', 'r', 'g'};


%% Define simulation directories and load data
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/csd-presentation');
mkdir(save_dir);

simDirs = ["simID_50", "simID_51", "simID_52"];
legendText = [
    "-30 km", "-20 km", "-15 km"
];
saveFilename = 'wTabmu_eq_lin_csd';
timestep = 2500;
lineMarkers = {'none', 'none', 'none', 'none'};
simColors = {'k', 'b', 'g', 'r'};
xLim = [-30, 0];
stLegendFontSize = 9;

% simDirs = ["simID_1", "simID_2", "simID_18"];
% legendText = [
%     "equilibrium",...
%     "linear ($\tau = 7$ days)",...
%     "csd (Hammer, Rutherford)"
% ];
% saveFilename = 'wTabmu_eq_lin_csd';
% timestep = 2800;
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'k', 'b', 'g', 'r'};
% xLim = [-35, 0];
% stLegendFontSize = 9;

% simDirs = ["simID_10", "simID_14", "simID_16", "simID_18"];
% legendText = [
%     "$\beta = \beta_{eq}$",...
%     "$I_0=10^{9}$, $U_0=10^{-10}$",...
%     "$I_0=10^{9}$, $U_0=10^{-9}$",...
%     "$I_0=10^{9}$, $U_0=10^{-8}$"];
% saveFilename = 'compare_csd.pdf';
% timestep = 3000;
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'k', 'b', 'g', 'r'};
% xLim = [-35, 0];
% stLegendFontSize = 12;

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
figX = 300 / inch;
figY = 0.618 * figX;

fig = figure('Name', 'all', 'Units','inches', 'Position',[5, 2, figX, figY]);
sgtitle(fig,  sprintf('Time = %0.1fh', FractureData{1}.time / 3600));
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
    alpha     = FractureData{i}.alpha;
    beta      = FractureData{i}.beta;
    xh20d      = FractureData{i}.xh2od;


    %% Width
    ax1 = subplot(2,2,1);
    X = xc(ext_elem);
    Y = Width(ext_elem);
    h_w = plot(ax1, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    xlabel(ax1, 'Depth (km)', FontSize=stFontSize);
    ylabel(ax1,"Width (m)", FontSize=stFontSize);
    grid on;
    set(ax1, xlim=xLim);
    title(ax1, "\bf (a)", FontSize=stFontSize);
    hold(ax1, 'on');


    %% Alpha and beta
    ax2 = subplot(2,2,2);
    X = xc(frac_elem);
    Alpha = alpha(frac_elem);
    h_alpha = plot(ax2, X, Alpha, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:20:length(X)); 
    hold(ax2, 'on');
    Beta = beta(frac_elem);
    h_beta = plot(ax2, X, Beta, Color=simColors{i}, LineStyle=lineStyle{2}, Marker=lineMarkers{i}, MarkerIndices=1:20:length(X));
    grid on;

    xlabel(ax2, 'Depth (km)', FontSize=stFontSize);
    ylabel(ax2,"Volume concentrations", FontSize=stFontSize);
    set(ax2, xlim=xLim);
    set(ax2, 'YTick',[0, 0.2, 0.4, 0.6, 0.8, 1]);
    title(ax2, "\bf (b)", FontSize=stFontSize);
    hold(ax2, 'on');

    %% Temperature
    ax3 = subplot(2,2,3);
    X = xc(frac_elem);
    Y = T(frac_elem);
    h_T = plot(ax3, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    grid on;
    xlabel(ax3, 'Depth (km)', FontSize=stFontSize);
    ylabel(ax3,"Temperature (C$^\circ$)", FontSize=stFontSize);
    set(ax3, xlim=xLim);
    title(ax3, "\bf (c)", FontSize=stFontSize);
    hold(ax3, 'on');

    %% Temperature
    ax4 = subplot(2,2,4);
    X = xc(frac_elem);
    Y = Viscosity(frac_elem);
    h_mu = semilogy(ax4, X, Y, Color=simColors{i}, LineStyle=lineStyle{1}, Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    grid on;
    xlabel(ax4, 'Depth (km)', FontSize=stFontSize);
    ylabel(ax4,"Viscosity (Pa$\cdot$s)", FontSize=stFontSize);
    set(ax4, xlim=xLim);
%     set(ax4, 'YTick',[1, 1e4, 1e8, 1e12, 1e16, 1e20],'YTickLabel',{'$1$', '$10^{4}$', '$10^{8}$', '$10^{12}$', '$10^{16}$', '$10^{20}$'});
    set(ax4, 'YTick',[1, 1e2, 1e4, 1e6, 1e8, 1e10],'YTickLabel',{'$1$', '$10^{2}$', '$10^{4}$', '$10^{6}$', '$10^{8}$', '$10^{10}$'});
%     set(ax4, 'YTick',[1, 1e1, 1e2, 1e3, 1e4, 1e5],'YTickLabel',{'$1$', '$10^{1}$', '$10^{2}$', '$10^{3}$', '$10^{4}$', '$10^{5}$'});
    title(ax4, "\bf (d)", FontSize=stFontSize);
    hold(ax4, 'on');
    
    hLine1 = [hLine1, h_w];
    if i == 1
        hLine2 = [h_alpha, h_beta];
    end
end
legend(ax4, hLine1, legendText, Location='northwest', FontSize=stLegendFontSize)
legend(ax2, hLine2, "bubble fraction, $\alpha$", "crystal content, $\beta$", Location='northwest', FontSize=stLegendFontSize)
% save_path = fullfile(save_dir, saveFilename);
% exportgraphics(fig, save_path + ".pdf");
% savefig(save_path);


function rho = rhoLith(x_km)
    alpha=0.001918;
    d=1.7247;
    rho1=2.7323;
    rho0=2.0929;
    rho = 1e3 * (rho1 - alpha * x_km - (rho1 - rho0) * exp(x_km / d));
end