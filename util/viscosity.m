clear all;
close all;
clc;
addpath(genpath(pwd));
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
set(gca, 'FontName', 'CMU Serif')

T0 = 850+273.15;
Lm = 350000;
GasDataPath      = 'StateEquation/propH2OCO2.mat'; % from VESICAL package
DissolvedGasPath = {
    'StateEquation/dpath0.1.mat',...
    'StateEquation/dpath0.4.mat',...
    'StateEquation/dpath0.6.mat',...
    'StateEquation/dpath0.8.mat'};   % state data for dissolved gas

States = cellfun(@(path) State(GasDataPath, path, T0, Lm), DissolvedGasPath, 'UniformOutput', false);

pressure = 1e6 * linspace(0.0, 1000.0, 100);
temperature = T0 * ones(1, length(pressure));
cellfun(@(state) state.updateDensityAndViscosity(pressure, temperature), States);

for i = 1:length(States)
    state = States{i};
    semilogy(1e-6 * pressure, state.muFrac, 'linewidth', 2)
    hold on;
end

xlabel('pressure, MPa', 'fontsize',14);
ylabel('viscosity, Pa$\cdot$s', 'fontsize',14);

legendText = cellfun(@(state) sprintf('H$_2$O saturation $=%0.1f \\%%$', 100*state.CH2O_sat), States, 'UniformOutput', false);
legend(legendText, 'Interpreter','latex', 'location','best', 'fontsize',14);

% set(gcf, 'InvertHardcopy', 'off')
saveas(gcf, [pwd '/images/viscosity_by_saturation.pdf'])
saveas(gcf, [pwd '/images/viscosity_by_saturation.png'])
% exportgraphics(gcf, [pwd '/images/viscosity_by_saturation.pdf'])