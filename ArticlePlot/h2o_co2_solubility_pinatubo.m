clear all;
close all;
clc;
addpath(genpath(pwd));

%% Redefine default fig parameteres to latex
stFontSize = 14;
textFontSize = 10;
set(0, 'defaultTextInterpreter','latex')
set(0, 'defaultTextFontName', 'CMU Serif')
set(0, 'defaultAxesFontName', 'CMU Serif')
set(0, 'defaultAxesFontSize', stFontSize)
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
set(0, 'defaultLineLineWidth', 2)
set(0, 'defaultLegendFontName', 'CMU Serif')
set(0, 'defaultLegendFontSize', stFontSize)
set(0, 'defaultLegendInterpreter', 'latex')

save_dir = fullfile(pwd, 'images/Article2022-review');
mkdir(save_dir);
isobars_data = importdata(fullfile(pwd, "ArticlePlot/pinatubo-isobars.csv"));
N = 50;
isobarsX = {};
isobarsY = {};
isobarsP = {200, 400, 600, 800, 900};
for i = 1:5
    r1 = (i-1)*N+1;
    r2 = i * N;
    x = isobars_data.data(r1:r2, 2);
    isobarsX{end+1} = isobars_data.data(r1:r2, 2)';
    isobarsY{end+1} = isobars_data.data(r1:r2, 3)';
end


isopleth_data = importdata(fullfile(pwd, "ArticlePlot/pinatubo-isopleth.csv"));
N = 50;
isoplethX = {};
isoplethY = {};
isoplethch2o = {0.2, 0.4, 0.6};
for i = 1:3
    r1 = (i-1)*N+1;
    r2 = i * N;
    x = isopleth_data.data(r1:r2, 2);
    isoplethX{end+1} = isopleth_data.data(r1:r2, 2)';
    isoplethY{end+1} = isopleth_data.data(r1:r2, 3)';
end

%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 160 / inch;
figY = 0.618 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);

for i = 1:5
    X = isobarsX{i};
    Y = isobarsY{i};
    plot(ax, X, Y, LineWidth=2, LineStyle='--', Color='black');
    hold on;
end
xlim([0 21])
ylim([0 0.6])
xlabel("$H_2O$ content, wt\%", FontSize=stFontSize)
ylabel("$CO_2$ content, wt\%", FontSize=stFontSize)


Pdpath=linspace(1, 900, 100);
load 'pdpath0.2.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h1 = plot(ax, X, Y, LineWidth=2, Color='b', Marker='o', MarkerIndices=1:5:length(X));

load 'pdpath0.4.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h2 = plot(ax, X, Y, LineWidth=2, Color='black', Marker='^', MarkerIndices=1:5:length(X));

load 'pdpath0.6.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h3 = plot(ax, X, Y, LineWidth=2, Color='r', Marker='s', MarkerIndices=1:5:length(X));

legendText = ["$c_{H_2O}\!=\!3.85$ wt\%", "$c_{H_2O}\!=\!6.23$ wt\%", "$c_{H_2O}\!=\!9.62$ wt\%"];
legend([h1, h2, h3], legendText, Location='best', FontSize=stFontSize)

xleft = [-1, 0.2, 1, 0, -0.5];
for i = 1:5
    X = isobarsX{i}+xleft(i);
    Y = isobarsY{i};
    xp = X(ceil(end/2));
    yp = Y(ceil(end/2));
    strText = num2str(isobarsP{i});
    text(ax, xp, yp, strText, FontSize=textFontSize, EdgeColor='k', BackgroundColor='white', HorizontalAlignment='center');
end
save_path = fullfile(save_dir, 'Solubility.pdf');
exportgraphics(fig, save_path);

% 
% for pp=20:100:1000
%     for j=1:101
%     hh(j)=h2o(Pd,Xd(j));
%     cc(j)=co2(Pd,Xd(j));
%     end
% plot(hh,cc,'LineWidth',2)
% % if(pp <1000) ,text(hh(end/2),cc(end/2),num2str(pp),'fontweight','bold','fontsize',14,'EdgeColor','k','BackgroundColor',[.7 .9 .7]);end
% 
% hold on
% end
% 
% for xx=0.:0.1:1
%     k=1;
%     for pp=50:100:1000
%         
%     hhh(k)=h2o(pp,xx);
%     ccc(k)=co2(pp,xx);
%     k=k+1;
%     end
%     
% plot(hhh(1:k-1),ccc(1:k-1),'--','LineWidth',1)
% 
% text(hhh(end),ccc(end)+0.015,num2str(xx),'fontweight','bold','fontsize',12)%,'EdgeColor','k','BackgroundColor',[.7 .9 .7]);
% hold on
% end
% 
% 
% xlabel('H_2O')
% ylabel('CO_2')
% set(gca,'Fontsize',16,'LineWidth',2,'fontweight','bold');
% 
% 
% ddd
% [pd,xd]=ndgrid(Pd,Xd);
% surf(pd,xd,h2o(pd,xd))
% shading interp
% xlabel('Pressure, MPa')
% ylabel('X_{H2O}')
% zlabel('H_2O,wt%')
% set(gca,'Fontsize',16,'LineWidth',2,'fontweight','bold');
% 
% 
% % surf(pp,xx,CO2);