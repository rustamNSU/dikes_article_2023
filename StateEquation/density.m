clear all
close all
clc
global mm0 mw0 mc0 h2o co2 p MH2O MCO2
MH2O=18.01528;
MCO2=44.01;
rhom0=2360.;
rhow0=852.3832924;
rhoc0=rhom0*0.9;
load('propH2OCO2.mat')
load('solub.mat')
h2o=griddedInterpolant(pp,xx,H2O);
co2=griddedInterpolant(pp,xx,CO2);
rho=griddedInterpolant(P,X,T,Rho);
xmol=0.3;
cw0=h2o(850,xmol)/100;
cc0=co2(850,xmol)/100;
ct0=cw0+cc0;
rho_md0=1/((1-ct0)/rhom0+cw0/rhow0+cc0/rhoc0);
mm0=rho_md0*(1-ct0);
mw0=rho_md0*cw0;
mc0=rho_md0*cc0;
wt_dry = 0.01; %H2O wt. % that defines "dry"
%[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
Porph = [70.61 0.28 17.43 1.44 0.01 0.33 1.95 5.06 2.79 0.04 wt_dry 0];
Composition = Porph; %Tuffen and Castro (2009) locality AO
T=850+273;

figure()
pause(1e-3)
px=850:-25:50;
den=px;
vis=px;
alph=px;
k=1;
for p=px
xW0=xmol*MH2O/(xmol*MH2O+(1-xmol)*MCO2);
options = optimset('Display','off');
xW=fsolve(@findXw,xW0,options);
MA=xW/MH2O+(1-xW)/MCO2;
xmol=xW/MH2O/MA;
cw=h2o(p,xmol)/100;
cc=co2(p,xmol)/100;
rho_g=rho(p,xW,T);
alpha=rhow0 * rhoc0 * ((cc + cw) * mm0 + (mc0 + mw0) * (-1 + cw + cc)) * rhom0 / (((((cc + cw) * mm0 + (mc0 + mw0) * (-1 + cw + cc)) * rhom0 + mm0 * rho_g * (-1 + cw + cc)) * rhow0 - cw * mm0 * rho_g * rhom0) * rhoc0 - cc * mm0 * rho_g * rhom0 * rhow0);
alph(k)=alpha;
rho_md=1/((1-cw-cc)/rhom0+cw/rhow0+cc/rhoc0);
den(k)=alpha*rho_g+(1-alpha)*rho_md;
vis(k) = Giordano_2008_visc(cw*100, T, Composition);
k=k+1;
pause(1e-8)
end

save('denvis.mat', 'px','den','vis')


yyaxis left
plot(px,den,'-','LineWidth',3);
xlabel('Pressure, MPa')
ylabel('Density, kg/m^3')
yyaxis right
semilogy(px,vis,'-', 'LineWidth', 3);
%ylabel('Bubble vol. fraction')
ylabel('Viscosity, Pa s')
xlim([0 850])
set(gca, 'LineWidth',1.5, 'FontSize', 16)
function res=findXw(xw)
global mm0 mw0 mc0 h2o co2 p MH2O MCO2
MA=xw/MH2O+(1-xw)/MCO2;
xmol=xw/MH2O/MA;
cw=h2o(p,xmol)/100;
cc=co2(p,xmol)/100;
XX=((-1 + cw + cc) * mw0 + cw * mm0) / ((-1 + cw + cc) * mw0 + (mc0 + mm0) * cw + (mc0 + mm0) * cc - mc0);
res=xw-XX;
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
    
%     for iz = 1:nx                     % step through each composition 
%     BT          = sum(bb.*bcf(iz,:));
%     CT          = sum(cc.*ccf(iz,:));
%     TG          = BT/(12-AT) + CT;
%     F           = BT/(TG*(1 - CT/TG)*(1 - CT/TG));
%     outvalues   =[outvalues ; iz AT BT CT TG F];
%     end
%     
%Calculate the coefficients using matrix algebra instead of a loop
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

