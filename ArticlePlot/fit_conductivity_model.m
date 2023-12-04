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

save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);

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
% Tr2 = Treservoir(Z, dTdz, omega, 950, D1, -32500);
p = reservoirPressure(Z);
rho = reservoirDensity(Z);
xh2o = 1.0;
hLines = [];
[Tl, Ts] = TlTs(p, xh2o);
hTs = plot(ax1, Ts, Z / 1e3, LineStyle='-', LineWidth=2, Color='b');   
hold on;
hTl = plot(ax1, Tl, Z / 1e3, LineStyle='-', LineWidth=2, Color='r');
hTr = plot(ax1, Tr, Z / 1e3, LineStyle='-', LineWidth=2, Color='k');
% hTr2 = plot(ax1, Tr2, Z / 1e3, LineStyle='-', LineWidth=2, Color='g');
xlabel("Temperature (C$^\circ$)");
ylabel("Depth (km)");
lgd1 = legend([hTr, hTl, hTs], ["reservoir", "liquidus", "solidus"], Location='southwest', FontSize=textFontSize, Interpreter='latex');
% title(lgd1, "$\mathbf{X_{H_2O}}$");
grid on;

ax2 = subplot(9,2,linspace(4,18,8));
t = readtable('ArticlePlot/wen_2015_granite_conductivity.csv')
Tf = table2array(t(:, 1));
kf = table2array(t(:, 2));

T = linspace(0, 1000, 100)
k0 = 1.9;
k = reservoirConductivity(T, k0);

plot(Tf, kf, 'o', Color='red', MarkerSize=5);
hold on;
plot(T, k, LineStyle='-', LineWidth=2, Color='k');
% save_path = fullfile(save_dir, 'reservoir_data.pdf');
% exportgraphics(fig, save_path, 'ContentType', 'vector');


function [Tl,Ts]=TlTs(PMPa,XH2Od)
    P=PMPa/1e6;
    dTliq=70;
    Tl=1060+0.18*P-10.2*sqrt(P.*XH2Od)-dTliq;
    Ts=Tl-(217+0.07*P)+dTliq;
end

function beq=betaeq(x)
    a = -2.836;  b = 22.12; c = -67.97;  d = 119.6;  f = -113.6; g = 45.78;
    beq=1./(1+exp(a+b*x+c*x.^2+d*x.^3+f*x.^4+g*x.^5));
end

function T=Treservoir(depth, dTdz, omega, Td, D1, chamberDepth)
    z = -depth;
    b=2e-3;
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

function k = reservoirConductivity(T, k0)
    b = 2e-3;
    k = k0 ./ (1.0 + b * T);
end