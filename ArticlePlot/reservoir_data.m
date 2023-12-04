clear all
close all
addpath(genpath(pwd));

% https://www.mathworks.com/matlabcentral/answers/597745-how-to-add-a-second-upper-x-axis-and-control-the-value-s-locations-and-text
%% Redefine default fig parameteres to latex
stFontSize = 12;
textFontSize = 12;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', textFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 250 / inch;
figY = 0.45 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);

ax1 = subplot(9,2,linspace(3,17,8));
Z = linspace(-40e3, 0.0, 100);
dTdz = 40e-3;
omega = 0.1;
Td = 850;
D1 = 10000;
Tr = Treservoir(Z, dTdz, omega, Td, D1, -32500);
p = reservoirPressure(Z);
rho = reservoirDensity(Z);

hTr = plot(ax1, Tr, Z / 1e3, LineStyle='-', LineWidth=2, Color='k');
xlabel("Temperature (C$^\circ$)");
ylabel("Depth (km)");
grid on;

ax2 = subplot(9,2,linspace(4,18,8));
xlabels{1} = 'Density (kg$/$m$^3$)';
xlabels{2} = 'Pressure (MPa)';
ylabels{1} = 'Depth (km)';
ylabels{2} = '';
[hrho, hp] = plotxx(ax2, rho, 1e-3*Z, 1e-6*p, 1e-3*Z, xlabels, ylabels);
lgd2 = legend(ax2, [hrho, hp], ["density", "lithostatic pressure"], Location='southwest', FontSize=textFontSize, Interpreter='latex');
ax2.YGrid = 'on';

function T=Treservoir(depth, dTdz, omega, Td, D1, chamberDepth)
    z = -depth;
    b=1.5e-3;
    H=-chamberDepth;  
    t2 = (D1 ^ 2);
    t4 = 0.1e1 / D1;
    t6 = exp(-t4 * H);
    t15 = log(0.1e1 / (Td * b + 1));
    t19 = exp(-t4 * z);
    t22 = z .^ 2;
    t24 = 2 * t2;
    t39 = H ^ 2;
    t46 = exp(0.1e1 ./ (0.2e1 * D1 * H + omega * t39 + 0.2e1 * t6 * t2 - t24) * ...
        (0.2e1 * t6 * dTdz * t2 * b * z - 0.2e1 * t19 * t2 * (H * b * dTdz + t15) + ...
        t15 * (-0.2e1 * D1 * z - t22 * omega + t24) + (H * omega * z + t24) .* (H - z) * b * dTdz));
    T = 0.1e1 / b * (t46 - 0.1e1);
end


% get density of reservoir (z is a depth in [m] with minus sign)
function rho = reservoirDensity(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    d = 1724.7;
    alpha = 0.001918;
    rho = rho1 - alpha * z - (rho1 - rho0) * exp(z / d);
end

% get lithostatic pressure of reservoir (z is a depth in [m] with minus sign)
function p = reservoirPressure(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    g = 9.81;
    d = 1724.7;
    alpha = 0.001918;
    p = g * (-rho1*z + 0.5*alpha*z.*z - d*(rho1-rho0)*(1.0-exp(z/d)));
end