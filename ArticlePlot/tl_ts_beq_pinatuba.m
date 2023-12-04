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

save_dir = fullfile(pwd, 'images/Article2022-review');
mkdir(save_dir);

%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 220 / inch;
figY = 0.618 * figX;
figY = 0.4 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax1 = subplot(1,2,1);

P = linspace(0,900) * 1e6;
xh2od = [0 0.25 0.5 0.75 1];
sio2 = 64.6;
% colors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE'};
% colors = {'#d7191c','#fdae61','#ffffbf','#abdda4','#2b83ba'};
hLines = [];
hDashedLines = [];
for i=1:5
    [TL, TS] = TlTS(P, xh2od(i), sio2);
    hDash = plot(ax1, TS-273.15, 1e-8*P, LineStyle='--', LineWidth=2, Color=colors{i});   
    hold on;
    h = plot(ax1, TL-273.15, 1e-8*P, LineStyle='-', LineWidth=2, Color=colors{i});
    hLines = [hLines h];
    hDashedLines = [hDashedLines hDash];
end
xlabel("Temperature (C$^\circ$)")
ylabel("Pressure (kbar)")
legendText = string(xh2od);
ylim([0,9]);
xlim([600, 1500]);
lgd1 = legend(hLines, legendText, Location='none', FontSize=textFontSize, Interpreter='latex', Position=[0.355 0.2 0.05 0.2]);
title(lgd1, "$\mathbf{X_{H_2O}}$");
grid on;

yyaxis(ax1,'right');
    z = [-5, -10, -15, -20, -25, -30, -35];
    p = 1e-8 * reservoirPressure(z*1e3);
    set(ax1, 'YColor','black');
    ylabel(ax1,"Depth (km)", FontSize=stFontSize);
    set(ax1, 'YTick',p/9, 'YTickLabel',-z);

ax2 = subplot(1,2,2);
fieldData = importdata('TX.csv');
Tred = fieldData.data(:,1);
X = 1 - fieldData.data(:,2);
Tr = linspace(0,1);

hField = plot(ax2, Tred, X, 'o', Color='red', MarkerSize=5);
hold on;
hFit = plot(ax2, Tr, betaeq(Tr), LineStyle='-', LineWidth=2, Color='black');
ylim([0 1])
xlim([0 1])
xlabel("Reduced temperature, $T'$")
ylabel('Crystal content, $\beta$')
legend([hField, hFit], ["experimental data", "curve fitting"], Location='best', FontSize=textFontSize)
grid on;
save_path = fullfile(save_dir, 'crystalization_by_p_and_t.pdf');
exportgraphics(fig, save_path, 'ContentType', 'vector');

% function [T_L,T_S]=TlTS(PMPa,XH2Od)
%     P=PMPa/100;
%     dTliq=70;
%     T_L=1060+18*P-102*sqrt(P.*XH2Od)-dTliq;
%     T_S=T_L-(217+7*P)+dTliq;
% end

% function beq=betaeq(x)
%     a = -2.836;  b = 22.12; c = -67.97;  d = 119.6;  f = -113.6; g = 45.78;
%     beq=1./(1+exp(a+b*x+c*x.^2+d*x.^3+f*x.^4+g*x.^5));
% end

function [TL, TS] = TlTS(P, XH2Od, SiO2)
    P = P / 1e8; % Pa to kbar
    aS=854.0896; bS=6; cS=224.0896; dS=80; eS=0.357; fS=6;
    aL=1205.7; bL=6; cL=285.7; dL=200; eL=0.7; fL=11;
    A=1287.6; B=-20.154;
    
    dTliq=A+B*SiO2;
    TS = 273.15 + aS+bS.*P-XH2Od.*(cS+fS.*P-dS./(P+eS));
    TL = 273.15 + aL+bL*P-XH2Od.*(cL+fL.*P-dL./(P+eL))+dTliq;
end

function beq = betaeq(x)
    aF=-4.974; bF=28.623; cF=-52.708; dF=34.816;
    beq = (1+exp(aF+bF.*(x)+cF.*(x).^2+dF.*(x).^3)).^(-1);
end

function p = reservoirPressure(z)
    rho0 = 2092.9;
    rho1 = 2732.3;
    g = 9.81;
    d = 1724.7;
    alpha = 0.001918;
    p = g * (-rho1*z + 0.5*alpha*z.*z - d*(rho1-rho0)*(1.0-exp(z/d)));
end