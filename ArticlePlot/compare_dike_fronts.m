clear all;
close all;
clc;
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
set(0, 'defaultLegendFontSize', stLegendFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

lineStyle = ["-", "-", "-", "-"];
lineMarkers = ["none", "none", "none", "none", "none", "none", "none", "none"];
xLim = [0.1, 100];
sim_root_dir = fullfile(pwd, 'Simulations');
generate_filename = @(i) sprintf('timestep_%05d.mat', i);
save_dir = fullfile(pwd, 'images/Article2022-review');
mkdir(save_dir);


% sim_dirs = {'simID_10011', 'simID_10021', 'simID_10031',...
%     'simID_10111', 'simID_10121', 'simID_10131',...
%     'simID_10211', 'simID_10221', 'simID_10231'};
% % legendText = {'simID 10011', 'simID 10021', 'simID 10031',...
% %               'simID 10111', 'simID 10121', 'simID 10131',...
% %               'simID 10211', 'simID 10221', 'simID 10231'};
% legendText = {
%     "k=0 [J/(m^2*K)], cH_2O=2.35%", "k=0 [J/(m^2*K)], cH_2O=6.24%", "k=0 [J/(m^2*K)], cH_2O=8.66%",...
%     "k=3 [J/(m^2*K)], cH_2O=2.35%", "k=3 [J/(m^2*K)], cH_2O=6.24%", "k=3 [J/(m^2*K)], cH_2O=8.66%",...
%     "k=5 [J/(m^2*K)], cH_2O=2.35%", "k=5 [J/(m^2*K)], cH_2O=6.24%", "k=5 [J/(m^2*K)], cH_2O=8.66%"
% };
% timesteps = {0:1:370, 0:1:114, 0:1:130,...
%              0:1:370, 0:1:370, 0:1:370,...
%              0:1:370, 0:1:228, 0:1:188};

% %% Heat transfer coefficient
% simDirs = ["simID_10021", "simID_10121", "simID_10221"];
% legendText = ["$k=0$", "$k=3$", "$k=5$"];
% saveFilename = 'compare_k_influence.pdf';
% timesteps = {0:1:114, 0:1:92, 0:1:81};
% 
% %% Water saturation
% simDirs = ["simID_10111", "simID_10121", "simID_10131"];
% legendText = ["$c_{H_2O}\!=\!2.35\%$", "$c_{H_2O}\!=\!6.24\%$", "$c_{H_2O}\!=\!8.66\%$"];
% saveFilename = 'compare_ch2o_influence.pdf';
% timesteps = {0:1:200, 0:1:92, 0:1:77};

%% Mass rate in chamber
% simDirs = ["simID_9121", "simID_10121", "simID_11121"];
% legendText = ["$Q_{ch}\!=\!1$", "$Q_{ch}\!=\!2$", "$Q_{ch}\!=\!4$"];
% saveFilename = 'compare_mass_rate.pdf';
% timesteps = {0:1:98, 0:1:92, 0:1:86};

% %% Different T0
% simDirs = ["simID_10121", "simID_20121"];
% legendText = ["$T_{ch}\!=\!850$ C$^\circ$", "$T_{ch}\!=\!900$ C$^\circ$"];
% saveFilename = 'compare_chamber_temperature.pdf';
% timesteps = {0:1:92, 0:1:92};
% simColors = {'black', 'r'};
% lineMarkers = {'^', 's'};

% %% Different k1c
% simDirs = ["simID_10121", "simID_31121"];
% legendText = ["$K_{1c}\!=\!1$ MPa$\cdot$m$^{1/2}$", "$K_{1c}\!=\!100$ MPa$\cdot$m$^{1/2}$"];
% saveFilename = 'compare_reservoir_toughness.pdf';
% timesteps = {0:1:92, 0:1:92};
% simColors = {'black', 'r'};
% lineMarkers = {'^', 's'};

%% Different E
% simDirs = ["simID_10121", "simID_32121"];
% legendText = ["$E\!=\!15$ GPa", "$E\!=\!45$ GPa"];
% saveFilename = 'compare_reservoir_young_front.pdf';
% timesteps = {0:1:92, 0:1:92};
% simColors = {'black', 'r'};
% lineMarkers = {'^', 's'};

% %% Compare tau
% simDirs = {'simID_40100', 'simID_40111', 'simID_40112'};
% legendText = ["equlibrium crystallization", "$\tau=100000$ s", "$\tau=604800$ s"];
% saveFilename = "images/noneqilibriumCryst/tau_variation_front";
% timesteps = {0:2:600, 0:2:600, 0:2:600};
% lineStyle = ["-", "-", "-", "-"];
% lineMarkers = ["none", "none", "none", "none", "none", "none", "none", "none"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Convergence Nx
% simDirs = ["simID_110002", "simID_110000", "simID_110001"];
% timesteps = {1:1:625, 1:1:625, 1:1:625};
% legendText = [
%     "$N_x=250$", ...
%     "$N_x=500$", ... 
%     "$N_x=1000$"];
% saveFilename = 'front_convergence_Nx.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r', 'g'};

% %% Conductivity
% simDirs = ["simID_110005", "simID_110000", "simID_110006"];
% timesteps = {1:1:625, 1:1:625, 1:1:625};
% legendText = [
%     "$k=0.5$ W/(m$\cdot$K)", ...
%     "$k=1$ W/(m$\cdot$K)", ... 
%     "$k=2$ W/(m$\cdot$K)"];
% saveFilename = 'front_conductivity.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r', 'g'};

% %% ch2o
% simDirs = ["simID_110008", "simID_110000", "simID_110009"];
% timesteps = {1:1:625, 1:1:625, 1:1:625};
% legendText = [
%     "$c_{H_2O}\!=\!2.35\%$", ...
%     "$c_{H_2O}\!=\!6.24\%$", ...
%     "$c_{H_2O}\!=\!8.66\%$"];
% saveFilename = 'front_ch2o.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r', 'g'};

% %% Injected mass
% simDirs = ["simID_110010", "simID_110000", "simID_110011"];
% timesteps = {1:1:625, 1:1:625, 1:1:625};
% legendText = [
%     "$M = 3\times10^7$ kg", ...
%     "$M = 5\times10^7$ kg", ...
%     "$M = 7\times10^7$ kg"];
% saveFilename = 'front_injected_mass.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r', 'g'};

%% Reservoir temperature
% simDirs = ["simID_140011", "simID_140001", "simID_140010"];
% timesteps = {1:3:4000, 1:2:1500, 1:2:1500};
% legendText = [
%     "$c_{H_2O}\!=\!2.35\%$", ...
%     "$c_{H_2O}\!=\!6.24\%$", ...
%     "$c_{H_2O}\!=\!8.66\%$"];
% saveFilename = 'front_chamber_temperature.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r'};

% simDirs = ["simID_140012", "simID_140001", "simID_140013"];
% timesteps = {1:2:1500, 1:2:1500, 1:2:1500};
% legendText = [
%     "$M = 5\times10^7$ kg", ...
%     "$M = 10\times10^7$ kg", ...
%     "$M = 15\times10^7$ kg"];
% saveFilename = 'front_chamber_temperature.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r'};
% xLim = [0.1, 100];

% simDirs = ["simID_140001", "simID_140014"];
% timesteps = {1:2:1500, 1:2:1500};
% legendText = ["$T_{ch}\!=\!850$ C$^\circ$", "$T_{ch}\!=\!900$ C$^\circ$"];
% saveFilename = 'front_chamber_temperature.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'k', 'r'};
% xLim = [0.1, 100];

% simDirs = ["simID_140015", "simID_140001", "simID_140016"];
% timesteps = {1:2:1500, 1:2:1500, 1:2:1500};
% legendText = ["$t_{inj}\!=\!20000$ s", "$t_{inj}\!=\!10000$ s", "$t_{inj}\!=\!5000$ s"];
% saveFilename = 'front_chamber_temperature.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'b', 'k', 'r'};
% xLim = [0.1, 100];


% simDirs = ["simID_142001", "simID_142010", "simID_142011", "simID_142012"];
% timesteps = {1:2:800, 1:2:800, 1:2:800, 1:2:800};
% legendText = [ 
%     "$K_{1c}=1$ MPa$\cdot$m$^{1/2}$", ...
%     "$K_{1c}=100$ MPa$\cdot$m$^{1/2}$", ...
%     "$K_{1c}=1000$ MPa$\cdot$m$^{1/2}$", ...
%     "$K_{1c}=1500$ MPa$\cdot$m$^{1/2}$"];
% saveFilename = 'front_toughness.pdf';
% lineMarkers = {'none', 'none', 'none', 'none'};
% simColors = {'k', 'b', 'g', 'r'};
% xLim = [0.1, 20];
% stLegendFontSize = 9;

% simDirs = ["simID_143017", "simID_143027", "simID_143037", "simID_143047", "simID_143057"];
% timesteps = {1:5:1000, 1:5:1000, 1:5:1000, 1:5:1000, 1:5:1000};
% legendText = [ 
%     "1", "2", "3", "4", "5"];
% saveFilename = 'compare_new_old.pdf';
% lineMarkers = {'none', 'none', 'none', 'none', 'none'};
% simColors = {'k', 'b', 'g', 'r', 'y'};
% xLim = [0.1, 400];
% stLegendFontSize = 9;


% simDirs = ["simID_143001", "simID_143102", "simID_143103"];
% timesteps = {1:4:2000, 1:4:2000, 1:4:2000};
% legendText = [
%     "$C_{SiO_2} = 64.6$ wt\%", ...
%     "$C_{SiO_2} = 70.6$ wt\%", ...
%     "$C_{SiO_2} = 70.6$ wt\%, $T_{ch}=758$ $^\circ$C"];
% saveFilename = 'compare_front_sio2.pdf';
% lineMarkers = {'none', 'none', 'none', 'none', 'none'};
% simColors = {'k', 'b', 'g', 'r', 'y'};
% xLim = [0.1, 20];
% stLegendFontSize = 9;

simDirs = ["simID_143105", "simID_143106"];
timesteps = {1:4:2000, 1:4:2000};
legendText = [
    "$C_{SiO_2} = 64.6$ wt\%", ...
    "$C_{SiO_2} = 70.6$ wt\%"];
saveFilename = 'compare_front_sio2_TL.pdf';
lineMarkers = {'none', 'none', 'none', 'none', 'none'};
simColors = {'k', 'b', 'g', 'r', 'y'};
xLim = [0.1, 20];
stLegendFontSize = 9;


sim_paths = fullfile(sim_root_dir, simDirs);
[front, time] = cellfun(@generate_front_depth_by_time, sim_paths, timesteps, 'UniformOutput', false);
[v, tv] = cellfun(@numeric_derivative, front, time, 'UniformOutput', false);

%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 150 / inch;
figY = 0.618 * figX;


fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);
yyaxis(ax,'left');
h_front = {};
for i = 1:length(simDirs)
    X = time{i} / 3600;
    Y = front{i} * 1e-3;
    h_1 = semilogx(ax, X, Y, Color=simColors{i}, LineStyle='-', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    h_front = [h_front, h_1];
    hold on;
end
ylim([-Inf, 0]);
set(ax, 'YColor','black');
ylabel(ax,"Depth, $l_t$ (km)", FontSize=stFontSize);

yyaxis(ax,'right');
for i = 1:length(simDirs)
    X = tv{i} / 3600;
%     X = 1:length(tv{i});
    Y = movmean(v{i},5);
    h_2 = loglog(ax, X, Y, Color=simColors{i}, LineStyle='--', Marker=lineMarkers{i}, MarkerIndices=1:10:length(X));
    hold on;
end
set(ax, 'xscale','log', 'yscale','log');
set(ax, 'YColor','black', 'YTick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],'YTickLabel',{'$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$'});
ylim([1e-5, 10]);
xlim(xLim);
ylabel(ax,"Velocity (m/s)", FontSize=stFontSize);
ax.XGrid = 'on';
xlabel(ax, 'Time (h)', FontSize=stFontSize)
set(ax, 'XTick',[1, 10, 100],'XTickLabel',{'$1$', '$10$', '$100$'});

% hStop = xline(ax, 10000/3600, Color='r', LineStyle='--', LineWidth=2);
legend(h_front, legendText, Location='west', FontSize=stLegendFontSize);

save_path = fullfile(save_dir, saveFilename);
exportgraphics(fig, save_path);