clear all
close all
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 12;
textFontSize = 9;
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
% figX = 183 / inch;
% figY = 0.618 * figX;
figX = 183 / inch;
figY = figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax1 = subplot(1,1,1);

Z = linspace(-37e3, 0.0, 100);
dTdz = 40e-3;
omega = 0.1;
Td = 850;
D1 = 10000;
Tr = Treservoir(Z, dTdz, omega, Td, D1, -32500);
P = Plith(Z);
xh2o = 1.0;
hLines = [];
[Tl, Ts] = TlTs(P, xh2o);
hTs = plot(ax1, Ts, Z / 1e3, LineStyle='-', LineWidth=2, Color='b');   
hold on;
hTl = plot(ax1, Tl, Z / 1e3, LineStyle='-', LineWidth=2, Color='r');
hTr = plot(ax1, Tr, Z / 1e3, LineStyle='-', LineWidth=2, Color='k');
xlabel("Temperature (C$^\circ$)")
ylabel("Depth (km)")

lgd1 = legend([hTr, hTl, hTs], ["reservoir", "liquidus", "solidus"], Location='southeast', FontSize=textFontSize, Interpreter='latex');
title(lgd1, "$\mathbf{X_{H_2O}}$")


% ax2 = subplot(1,2,2);
% fieldData = importdata('TX.csv');
% Tred = fieldData.data(:,1);
% X = 1 - fieldData.data(:,2);
% Tr = linspace(-0.2,1.2);

% hField = plot(ax2, Tred, X, 'o', Color='red', MarkerSize=5);
% hold on;
% hFit = plot(ax2, Tr, betaeq(Tr), LineStyle='-', LineWidth=2, Color='black');
% ylim([0 1])
% xlim([-0.2 1.2])
% xlabel("Reduced temperature, $T'$")
% ylabel('Crystal content, $\beta$')
% legend([hField, hFit], ["experimental data", "curve fitting"], Location='best', FontSize=textFontSize)

% save_path = fullfile(save_dir, 'crystalization_by_p_and_t.pdf');
% exportgraphics(fig, save_path, 'ContentType', 'vector');
% 
function [Tl,Ts]=TlTs(PMPa,XH2Od)
    P=PMPa/1e6;
    dTliq=70;
    Tl=1060+0.18*P-10.2*sqrt(P.*XH2Od)-dTliq;
    Ts=Tl-(217+0.07*P)+dTliq;
end
% function [T_L,T_S]=TlTs(PMPa,XH2Od)
%     P=PMPa/100;
%     dTliq=70;
%     T_L=1060+18*P-102*sqrt(P.*XH2Od)-dTliq;
%     T_S=T_L-(217+7*P)+dTliq;
% end

function beq=betaeq(x)
    a = -2.836;  b = 22.12; c = -67.97;  d = 119.6;  f = -113.6; g = 45.78;
    beq=1./(1+exp(a+b*x+c*x.^2+d*x.^3+f*x.^4+g*x.^5));
end

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

function Pl = Plith(depth)
    Z = 1e-3 * (0 - depth)
    alpha=0.001918;
    d=1.7247;
    rho1=2.7323;
    rho0=2.0929;
    g=9.81;
    Pl = (-0.2e1 * d * exp(-Z / d) * (-rho1 + rho0) + (0.2e1 * rho0 - 0.2e1 * rho1) * d ...
        + Z .* (Z * alpha + 0.2e1 * rho1)) * g / 0.2e1;
    Pl = 1e6 * Pl;
end