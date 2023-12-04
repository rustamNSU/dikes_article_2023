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

save_dir = fullfile(pwd, 'images/Article2022');
mkdir(save_dir);

res = importdata(fullfile(pwd, 'StateEquation', 'largetable.csv'));
MH2O = 18.01528;
MCO2 = 44.01;

P=res(:,1);
XH2O=res(:,3);
XX=XH2O*MH2O./(XH2O*MH2O+(1-XH2O)*MCO2);
DH2O=res(:,4);
DCO2=res(:,5); 

pp=reshape(P,21,11)/10;
xx=reshape(XX,21,11);
H2O=reshape(DH2O,21,11);
CO2=reshape(DCO2,21,11);
h2o=griddedInterpolant(pp,xx,H2O);
co2=griddedInterpolant(pp,xx,CO2);


%% Create figure 89mm for one column, 183 mm for two columns
inch = 25.4; %mm
figX = 160 / inch;
figY = 0.618 * figX;
fig = figure('Name', 'all', 'Units','inches', 'Position',[3, 2, figX, figY]);
ax = subplot(1,1,1);

Pd=50:100:950;
Xd=0:0.01:1;
[pg,xg]=ndgrid(Pd,Xd);
hhp=h2o(pg,xg);
ccp=co2(pg,xg);  
plot(ax, hhp',ccp', LineWidth=2, LineStyle='--', Color='black');
hold on
xlim([0 21])
ylim([0 0.48])
xlabel("$H_2O$ content, wt\%", FontSize=stFontSize)
ylabel("$CO_2$ content, wt\%", FontSize=stFontSize)


Pdpath=linspace(1, 950, 100);
load 'dpath0.1.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h1 = plot(ax, X, Y, LineWidth=2, Color='b', Marker='o', MarkerIndices=1:5:length(X));

load 'dpath0.4.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h2 = plot(ax, X, Y, LineWidth=2, Color='black', Marker='^', MarkerIndices=1:5:length(X));

load 'dpath0.6.mat'
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
X = h2o(Pdpath)*100;
Y = co2(Pdpath)*100;
h3 = plot(ax, X, Y, LineWidth=2, Color='r', Marker='s', MarkerIndices=1:5:length(X));

legendText = ["$c_{H_2O}\!=\!2.35$ wt\%", "$c_{H_2O}\!=\!6.24$ wt\%", "$c_{H_2O}\!=\!8.66$ wt\%"];
legend([h1, h2, h3], legendText, Location='best', FontSize=stFontSize)

for i=1:numel(Pd)
    xp = hhp(i,ceil(end/2));
    yp = ccp(i,ceil(end/2));
    strText = num2str(Pd(i));
    text(ax, xp, yp, num2str(Pd(i)), FontSize=textFontSize, EdgeColor='k', BackgroundColor='white', HorizontalAlignment='center');
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