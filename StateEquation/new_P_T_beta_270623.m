clear all
close all
figure('Units','normalized','Color','white','Position',[.10,.10,.8,.6])
ttt=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
global Pres Temp rho Tk cv uGas co2 h2o QH2O QCO2 Qc rho_c Composition beta rho_t np nt CH2O  Xf rhoGas rhoCrystal alphaFig betaFig

load('propH2OCO2.mat')
load('dpath0.6.mat')
Tk=0;
h2o=griddedInterpolant(PMPa,H2O_liq);
co2=griddedInterpolant(PMPa,CO2_liq);
Xf =griddedInterpolant(PMPa,XH2O_fl);
rho=griddedInterpolant(P,X,T-Tk,Rho);
cv=griddedInterpolant(P,X,T-Tk,Cv);
uGas = griddedInterpolant(P,X,T-Tk,U);

Composition = [70.61 0.28 17.43 1.44 0.01 0.33 1.95 5.06 2.79 0.04 0.01 0];

%%parameters
np=150;                                % #points in pressure array
nt=6;
Temp0=850;                               %temperature in C             range [800  900]
Psat=PMPa(end);                              %Saturatiom pressure in MPa
CH2O_sat=H2O_liq(end);                          %saturation magma water content range [0.03-0.1]
CCO2_sat=CO2_liq(end);                %saturation magma CO2 content
XH2Od=CH2O_sat/(CH2O_sat+CCO2_sat);
beta0=betaeq(Psat,XH2Od,Temp0);          %crystal content at saturation
rho_c=2700.;                            %density of the crystal phase
rho_m=Rhom(CH2O_sat,CCO2_sat);          %density of the melt phase at saturation

Qc  =rho_m*(1-beta0)*(1-CCO2_sat-CH2O_sat)+rho_c*beta0;   %density of condence phase at saturation
QH2O=rho_m*(1-beta0)*CH2O_sat;                            %density of dissolved water
QCO2=rho_m*(1-beta0)*CCO2_sat;                            %density of dissolved CO2

Presure=linspace(10,Psat,np);
Temperture(1:nt)=linspace(700+273.15,950+273.15,nt);
[Pres,Temp]=ndgrid(Presure,Temperture);
[rho_t,U_t,alpha,beta,CH2O,XH2O]=den(Pres,Temp);
Tu=scatteredInterpolant(Pres(:),U_t(:),Temp(:));


nexttile(ttt)
plot(Pres,alpha,'linewidth',2);
ylabel('Bubble fraction')
xlabel('Pressure, MPa')
set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)

nexttile(ttt)
set(gca,'visible','off')
% plot(Pres,XH2O,'linewidth',2);
% set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
% ylabel('XH2O')
% xlabel('Pressure, MPa')

nexttile(ttt)
plot(Pres,beta,'linewidth',2);
ylabel('Crystal content')
xlabel('Pressure, MPa')

set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)

nexttile(ttt)
plot(Pres,rho_t,'linewidth',2);
set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
ylabel('density kg/m^3')
xlabel('Pressure, MPa')


nexttile(ttt)
mu=CH2O;
for iP=1:np
    for iT=1:nt
        muu= Giordano_2008_visc(CH2O(iP,iT)'*100, Temperture(iT)+Tk, Composition)'.*theta(beta(iP,iT));
        mu(iP,iT)=muu;
    end
end
semilogy(Pres,mu,'linewidth',2);
ylabel('Viscosity, Pa\bullets')
xlabel('Pressure, MPa')
set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
hold on

nexttile(ttt)

plot(Pres,U_t/1e3./rho_t,'-','linewidth',2);
set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
ylabel('int. energy density, kJ/kg')
xlabel('Pressure, MPa')
legendStr = string(squeeze(Temperture(:)-273.15));
hl=legend(legendStr);
title(hl,'Temperature, ^oC');
hold on

% nexttile(ttt)
% plot(Pres,rhoGas,'linewidth',2);
% set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
% ylabel('gas density relative kg/m^3')
% xlabel('Pressure, MPa')
% 
% nexttile(ttt)
% plot(Pres,rhoCrystal,'linewidth',2);
% set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
% ylabel('crystal density relative kg/m^3')
% xlabel('Pressure, MPa')

% nexttile(ttt)
% plot(Pres,alphaFig,'linewidth',2);
% set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
% ylabel('alpha')
% xlabel('Pressure, MPa')
% 
% nexttile(ttt)
% plot(Pres,betaFig,'linewidth',2);
% set(gca,'FontSize',12, 'Fontweight', 'bold','linewidth',1.5)
% ylabel('beta')
% xlabel('Pressure, MPa')


function rho_m=Rhom(CH20,CCO2)
rhom0=2360.;
rhow0=852.;
rhoc0=rhom0*0.9;
rho_m=1./((1-CCO2-CH20)/rhom0+CH20/rhow0+CCO2/rhoc0);
end



function v = Giordano_2008_visc(H2Ot, T, Composition)
%Get the Giordano2008 VFT parameters from their script
%Giordano2008_Model.mat
outvalues=Giordano2008_Model(H2Ot, Composition);
At = outvalues(:,1);
Bt = outvalues(:,2);
Ct = outvalues(:,3);
v = 10.^(At+(Bt./(T-Ct)));
end

function outvalues=Giordano2008_Model(H2Ot, Composition)
%This is code is from Girodano 2008 (see citation below) modified
%for use by the Coumans et al., 2020 bubble growth model

% SCRIPT grdmodel08 Nov 2007
%
% MATLAB script to compute silicate melt viscosity.
%
% Citation: Giordano D, Russell JK, & Dingwell DB (2008)
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
%INPUT: Chemical compositions of silicate melts (Wt. % oxides) as:
% SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% One line for each melt composition
% _________________________________________________________________

% ________________________________________________________________
% OUTPUT: VFT Model parameters A, B and C for each melt composition
%         to model temperature dependence: log n = A + B/[T(K) - C]
% VFT_model.out contains: No. , A_value , B_value , C_value , Tg , F
% VFT_curves.out contains: Temperature Array in Kelvins &
%                          Log(n) values at T(K) values (1 line per melt)
% ________________________________________________________________

% VFT Multicomponent-Model Coefficients
% -----------------------------------------------------------------
AT      = -4.55;
bb  = [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  = [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

row = length(H2Ot);
Comp_rep = repmat(Composition,row,1);
Comp_rep(:,11) = H2Ot';

wtm=Comp_rep;

[nx,nc]=size(wtm);

bb_rep = repmat(bb,nx,1);
cc_rep = repmat(cc,nx,1);
AT_rep = repmat(AT,nx,1);

% Function molefrac_grd: converts wt % oxide bais to mole % oxide basis
%[ncomps xmf_t] = molefrac_grd(wtm);
[ncomps, xmf_t] = molepct_grd(wtm);

%outvalues=[];
% Load composition-basis matrix for multiplication against model-coefficients
% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
siti    =   xmf_t(:,1) + xmf_t(:,2);
tial    =   xmf_t(:,2)+xmf_t(:,3);
fmm     =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,6);
nak     =   xmf_t(:,8) + xmf_t(:,9);
b1  =   siti;
b2  =   xmf_t(:,3);
b3  =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,10);
b4  =   xmf_t(:,6);
b5  =   xmf_t(:,7);
b6  =   xmf_t(:,8) + xmf_t(:,11) + xmf_t(:,12);
b7  =   xmf_t(:,11) + xmf_t(:,12) + log(1+xmf_t(:,11));
b12 =   siti.*fmm;
b13 =   (siti + xmf_t(:,3) + xmf_t(:,10)).*( nak + xmf_t(:,11) );
b14 =   xmf_t(:,3).*nak;

c1      =   xmf_t(:,1);
c2      =   tial;
c3      =   fmm;
c4      =   xmf_t(:,7);
c5      =   nak;
c6      =   log(1+xmf_t(:,11) + xmf_t(:,12));
c11     =   xmf_t(:,3) + fmm + xmf_t(:,7) - xmf_t(:,10);
c11     =   c11.*(nak + xmf_t(:,11) + xmf_t(:,12));
bcf      =   [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf      =   [c1 c2 c3 c4 c5 c6 c11];

BT          = sum(bb_rep.*bcf(:,:),2);
CT          = sum(cc_rep.*ccf(:,:),2);
TG          = BT./(12-AT_rep) + CT;
F           = BT./(TG.*(1 - CT./TG).*(1 - CT./TG));

outvalues   =[AT_rep BT CT TG F];
end

function [nox xmf] = molepct_grd(wtm)
[nr nc] =   size(wtm);
% [n x]=MOLEPCT(X) Calculates mole percent oxides from wt %
% 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1
% Output: mole fractions of equivalent
mw=[60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];
mp=[];
xmf=[];

%Replicate the molecular weight row by the number of spatial nodes used
%(number of data rows)
mw_rep = repmat(mw,nr,1);

%Perform the calculation using matrix algebra
wtn = [wtm(:,1:10).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12)) wtm(:,11) 0.5.*wtm(:,12).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12))];
mp=wtn./mw_rep;
mpv= 100*(mp./sum(mp,2));
xmf=mpv;

[nr, nc]=size(xmf);
nox =nc;
end

function Theta=theta(bx)
be0=0.45;
teta1=10d0.^(1.884981564.*(bx-be0)+4.546330438.*(bx-be0).*(bx-be0));
teta1(bx<0.25)=0.63;
aa1= 0.999916;
aa2= 1.485884;
aa3= 3.98937;
aa4= 16.9386;
Fbx=erf(0.88622796./aa1.*(aa2*bx).*(1d0+(aa2*bx).^aa3));
teta2 = (1e0+(aa2*bx).^aa4).*(1d0-aa1.*Fbx).^(-2.5d0/aa2);
Theta=max(1.,teta1.*teta2);
end

function [rho_t,U_t,alpha,beta,CH2O,XH2O]=den(P,T)
global rho cv uGas co2 h2o Xf rho_c QH2O QCO2 Qc rhoGas rhoCrystal alphaFig betaFig Tk
CCO2=co2(P);
CH2O=h2o(P);
XH2O=Xf(P);
XH2Od=CH2O./(CH2O+CCO2);
beta=betaeq(P,XH2Od,T-273.15);
rho_m=Rhom(CH2O,CCO2);
rho_g=rho(P,XH2O,T);
% dt = 1;
U_g=uGas(P,XH2O,T);
% U_g1=uGas(P,XH2O,T-dt);
% Cv_g = cv(P,XH2O,T);
Cv_m=1200;
Cv_s=1200;

t1 = CCO2 + CH2O - 1;
t9 = rho_m .* (beta - 1) .* (QCO2 .* t1 + QH2O .* t1 + Qc .* (CCO2 + CH2O));
t12 = beta .* rho_c .* (QCO2 + QH2O);
alpha = max(0,1 ./ (Qc .* rho_g + t12 + t9) .* (t9 + t12));

rho_t=((beta - 1) .* rho_m - rho_c .* beta + rho_g) .* alpha + (1 - beta) .* rho_m + rho_c .* beta;

rhoGas = alpha .* rho_g;
rhoCrystal = rho_c.*(1-alpha).*beta;
rhoMelt = rho_m.*(1-alpha).*(1-beta);
alphaFig = alpha;
betaFig = beta;
U_t=((rhoCrystal .* Cv_s + rhoMelt .* Cv_m) .* (T+Tk) + rhoGas .* U_g)./ rho_t;
end

% function beq=betaeq(Pm,XH2Od,T)
% P=Pm/100;
% dTliq=70;
% T_L=1060+18*P-102*sqrt(P.*XH2Od)-dTliq;
% T_S=T_L-(217+7*P)+dTliq;
% x=(T-T_S)./(T_L-T_S);
% a = -2.836;  b = 22.12; c = -67.97;  d = 119.6;  f = -113.6; g = 45.78;
% beq=1./(1+exp(a+b*x+c*x.^2+d*x.^3+f*x.^4+g*x.^5));
% 
% 
% % beq=1-(2.92*Theta-4.76*Theta.^2+2.84*Theta.^3);
% % beq(T<T_S)=1;
% % beq(T>T_L)=0;
% % 
% end

function F=betaeq(Pm,XH2Od,Ti )
    aS=854.0896; bS=6; cS=224.0896; dS=80; eS=0.357; fS=6;
    aL=1205.7; bL=6; cL=285.7; dL=200; eL=0.7; fL=11; SiO2=70.2;
    aF=-4.974; bF=28.623; cF=-52.708; dF=34.816; A=1287.6;B=-20.154;
    P=Pm/100; dTliq=A+B*SiO2;
    TS = aS+bS.*P-XH2Od.*(cS+fS.*P-dS./(P+eS));
    TL = aL+bL*P-XH2Od.*(cL+fL.*P-dL./(P+eL))+dTliq;
    T = (Ti-TS)./(TL-TS);
    F = (1+exp(aF+bF.*(T)+cF.*(T).^2+dF.*(T).^3)).^(-1);
end