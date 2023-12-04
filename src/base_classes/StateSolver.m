classdef StateSolver < matlab.mixin.Copyable
    properties (SetAccess=public)
        pathToPropGas          % path to file where saved gas density by P and T
        pathToPropDissolvedGas % path to file where saved dissolved gas density by P and T
        Pvesical    % Pressure from VESIcal
        Tvesical    % Temperature from VESIcal
        T1d         % Temperature range
        Xvesical    % Mass fracton of water from VESIcal
        Rhog        % gas density by P, T, Xf from VESIcal
        Ugas
        h2o
        co2
        Xf
        oxides_composition
        
        Cv_c  % crystal specific heat (for unit mass)
        Cv_m  % melt specific heat
        Lc    % latent heat due crystallization
        tau   % [s], relaxation time in the crystallization law
        

        TempInChamber   % initial temperature (in chamber) [K]
        PInChamber
        uInChamber      % internal energy in chamber [J/kg]
        rhoInChamber    % density in chamber [kg/m^3]
        rhoCrystalInChamber
        rhom0   % melt phase density
        rhow0   % dissolved water density
        rhoCO20 % dissolved CO2 density
        rhoc0   % density of the crystal phase
        beta0
        Qc    % density of condence phase at saturation
        QH2O  % density of dissolved water
        QCO2  % density of dissolved CO2
    end

    methods (Access = private)
        function rho_m = Rhom(obj, cH2O, cCO2)
            c1 = 1-cCO2-cH2O;
            c2 = cH2O;
            c3 = cCO2;
            rho1 = obj.rhom0;
            rho2 = obj.rhow0;
            rho3 = obj.rhoCO20;
            rho_m = 1 ./ (c1 / rho1 + c2 / rho2 + c3 / rho3);
        end
    end

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = StateSolver(pathToPropGas, pathToPropDissolvedGas, magmaInput, T0)
            if nargin > 0
                obj.pathToPropGas = pathToPropGas;
                obj.pathToPropDissolvedGas = pathToPropDissolvedGas;
                load(pathToPropGas, 'P', 'Rho', 'T', 'X', 'U');
                load(pathToPropDissolvedGas, 'PMPa', 'H2O_liq', 'CO2_liq', 'XH2O_fl');
                obj.Lc = magmaInput.Lc;
                obj.Cv_c = magmaInput.Cv_c;
                obj.Cv_m = magmaInput.Cv_m;
                obj.tau = magmaInput.tau;
                
                obj.Pvesical = P;
                obj.Tvesical = T;
                obj.Xvesical = X;
                obj.Rhog = griddedInterpolant(P, X, T, Rho);
                obj.Ugas = griddedInterpolant(P, X, T, U);
                obj.T1d = squeeze(T(1, 1, :))';

                obj.h2o = griddedInterpolant(PMPa,H2O_liq);
                obj.co2 = griddedInterpolant(PMPa,CO2_liq);
                obj.Xf  = griddedInterpolant(PMPa,XH2O_fl);
                obj.oxides_composition = magmaInput.oxides_composition;


                obj.TempInChamber = T0;        % initial temperature (in chamber)
                obj.rhom0 = 2360.;     % melt phase density
                obj.rhow0 = 852.;      % dissolved water density
                obj.rhoCO20 = obj.rhom0*0.9; % dissolved CO2 density
                obj.rhoc0 = 2700.;     % density of the crystal phase

                Psat     = PMPa(end);    %Saturatiom pressure in MPa
                CH2O_sat = H2O_liq(end); %saturation magma water content range [0.03-0.1]
                CCO2_sat = CO2_liq(end); %saturation magma CO2 content
                XH2Od    = CH2O_sat / (CH2O_sat + CCO2_sat);

                SiO2 = obj.oxides_composition(1);
                [obj.beta0, TL, TS] = obj.betaeq(Psat * 1e6, XH2Od, T0, SiO2); % crystal content at saturation
                rho_m = obj.Rhom(CH2O_sat,CCO2_sat);    % density of the melt phase at saturation

                obj.Qc   = rho_m*(1-obj.beta0)*(1-CCO2_sat-CH2O_sat)+obj.rhoc0*obj.beta0; % density of condence phase at saturation
                obj.QH2O = rho_m*(1-obj.beta0)*CH2O_sat; % density of dissolved water
                obj.QCO2 = rho_m*(1-obj.beta0)*CCO2_sat; % density of dissolved CO2
            end            
        end%constructor
        
        function updateBetaChamber(obj, beta0)
            obj.beta0 = beta0;
            load(obj.pathToPropDissolvedGas, 'H2O_liq', 'CO2_liq');
            CH2O_sat = H2O_liq(end); %saturation magma water content range [0.03-0.1]
            CCO2_sat = CO2_liq(end); %saturation magma CO2 content
            rho_m = obj.Rhom(CH2O_sat,CCO2_sat);
            obj.Qc   = rho_m*(1-obj.beta0)*(1-CCO2_sat-CH2O_sat)+obj.rhoc0*obj.beta0; % density of condence phase at saturation
            obj.QH2O = rho_m*(1-obj.beta0)*CH2O_sat; % density of dissolved water
            obj.QCO2 = rho_m*(1-obj.beta0)*CCO2_sat; % density of dissolved CO2
        end

        function updateEqDensityAndViscosity(obj, newFrac, settings)
            fracElements = newFrac.get_fracture_elements();
            P = newFrac.pressure(fracElements); % convert Pa to MPa
            PMPa = P / 1e6;
            T = newFrac.temperature(fracElements);

            cCO2 = obj.co2(PMPa);
            cH2O = obj.h2o(PMPa);
            XH2O = obj.Xf(PMPa);
            XH2Od= cH2O ./ (cH2O + cCO2);

            newFrac.cco2(fracElements)  = cCO2;
            newFrac.ch2o(fracElements)  = cH2O;
            newFrac.xh2o(fracElements)  = XH2O;
            newFrac.xh2od(fracElements) = XH2Od;
            
            SiO2 = obj.oxides_composition(1);
            [betaeq, TL, TS]  = obj.betaeq(P, XH2Od, T, SiO2);
            newFrac.betaeq(fracElements) = betaeq; 
            beta = betaeq;

            rho_m = obj.Rhom(cH2O, cCO2);    % unit density of melt phase
            rho_g = obj.Rhog(PMPa, XH2O, T);    % unit density of gas
            
            t1 = cCO2 + cH2O - 1;
            t9 = rho_m .* (beta - 1) .* (obj.QCO2 .* t1 + obj.QH2O .* t1 + obj.Qc .* (cCO2 + cH2O));
            t12 = beta .* obj.rhoc0 .* (obj.QCO2 + obj.QH2O);
            alpha = max(0,1 ./ (obj.Qc .* rho_g + t12 + t9) .* (t9 + t12));

            newFrac.alpha(fracElements) = alpha;
            newFrac.beta(fracElements) = beta;
            newFrac.rhoc(fracElements) = (1.0 - alpha) .* beta .* obj.rhoc0;
            newFrac.rhog(fracElements) = alpha .* rho_g;
            newFrac.rhom(fracElements) = (1.0 - alpha) .* (1.0 - beta) .* rho_m;
            newFrac.rho(fracElements) = newFrac.rhog(fracElements) + newFrac.rhom(fracElements) + newFrac.rhoc(fracElements);
            
            newFrac.Mc(fracElements) = newFrac.rhoc(fracElements) ./ newFrac.rho(fracElements);
            newFrac.Mg(fracElements) = newFrac.rhog(fracElements) ./ newFrac.rho(fracElements);
            newFrac.Mm(fracElements) = newFrac.rhom(fracElements) ./ newFrac.rho(fracElements);

            obj.updateViscosity(newFrac);
        end
        
        % Calculate magma density by pressure and temperature
        % P and T can be number or array and has same length
        % P in [Pa]
        % T in [C]
        function updateDensity(obj, newFrac, oldFrac, startTime, endTime, mQ, settings)
            fracElements = newFrac.get_fracture_elements();
            P = newFrac.pressure(fracElements); % convert Pa to MPa
            PMPa = P / 1e6;
            PMPa(PMPa < 50) = 50;
            T = newFrac.temperature(fracElements);

            cCO2 = obj.co2(PMPa);
            cH2O = obj.h2o(PMPa);
            XH2O = obj.Xf(PMPa);
            XH2Od= cH2O ./ (cH2O + cCO2);

            newFrac.cco2(fracElements)  = cCO2;
            newFrac.ch2o(fracElements)  = cH2O;
            newFrac.xh2o(fracElements)  = XH2O;
            newFrac.xh2od(fracElements) = XH2Od;

            SiO2 = obj.oxides_composition(1);
            [betaeq, TL, TS] = obj.betaeq(P, XH2Od, T, SiO2);
            newFrac.TL(fracElements) = TL;
            newFrac.TS(fracElements) = TS;
            
            newFrac.betaeq(fracElements) = betaeq;
            beta = betaeq;
            rho_m = obj.Rhom(cH2O, cCO2);    % unit density of melt phase
            rho_g = obj.Rhog(PMPa, XH2O, T); % unit density of gas
            
            t1 = cCO2 + cH2O - 1;
            t9 = rho_m .* (beta - 1) .* (obj.QCO2 .* t1 + obj.QH2O .* t1 + obj.Qc .* (cCO2 + cH2O));
            t12 = beta .* obj.rhoc0 .* (obj.QCO2 + obj.QH2O);
            alpha = max(0,1 ./ (obj.Qc .* rho_g + t12 + t9) .* (t9 + t12));

            newFrac.alpha(fracElements) = alpha;
            newFrac.beta(fracElements) = beta;
            newFrac.rhoc(fracElements) = (1.0 - alpha) .* beta .* obj.rhoc0;
            newFrac.rhog(fracElements) = alpha .* rho_g;
            newFrac.rhom(fracElements) = (1.0 - alpha) .* (1.0 - beta) .* rho_m;
            newFrac.rho(fracElements) = newFrac.rhog(fracElements) + newFrac.rhom(fracElements) + newFrac.rhoc(fracElements);
            
            newFrac.Mc(fracElements) = newFrac.rhoc(fracElements) ./ newFrac.rho(fracElements);
            newFrac.Mg(fracElements) = newFrac.rhog(fracElements) ./ newFrac.rho(fracElements);
            newFrac.Mm(fracElements) = newFrac.rhom(fracElements) ./ newFrac.rho(fracElements);
        end

        % Calculate magma viscosity by pressure and temperature
        % T, cH2O, beta is number only !!!
        function updateViscosity(obj, newFrac)
            fracElements = newFrac.get_fracture_elements();
            for i=1:length(fracElements)
                element = fracElements(i);
                t = newFrac.temperature(element);
                ch2o = newFrac.ch2o(element);
                b = newFrac.beta(element);
                a = newFrac.alpha(element);
                newFrac.mu(element) = obj.Giordano_2008_visc(ch2o * 100, t, obj.oxides_composition) * obj.theta(b) * obj.bubble_viscosity(a);
            end
        end

        function updateDensityAndViscosity(obj, newFrac, oldFrac, startTime, endTime, mQ, settings)
            obj.updateDensity(newFrac, oldFrac, startTime, endTime, mQ, settings);
            obj.updateViscosity(newFrac);
        end
        
        
        function T = findTemperature(obj, frac)
            fracElements = frac.get_fracture_elements();
            T = zeros(size(fracElements));
            InternalEnergy = frac.u(fracElements);
            P = frac.pressure(fracElements) / 1e6;
            for i = 1:length(fracElements)
                element = fracElements(i);
                p = P(i) * ones(size(obj.T1d));
                x = frac.xh2o(element) * ones(size(obj.T1d));
                ug = squeeze(obj.Ugas(p, x, obj.T1d));
                u = frac.Mc(element) * obj.Cv_c * obj.T1d +...
                    frac.Mm(element) * obj.Cv_m * obj.T1d +...
                    frac.Mg(element) * ug;
                try
                    Tinterp = griddedInterpolant(u, obj.T1d);
                    T(i) = Tinterp(InternalEnergy(i));
                    T(i) = max(T(i), 600 + 273.15);
                    T(i) = min(T(i), 900 + 273.15);
                catch
                    a = 10;
                end
            end
            frac.temperature(fracElements) = T;
        end            
        

        function set_u0(obj, Pressure, Temperature)
            obj.PInChamber = Pressure;
            P = Pressure; % convert Pa to MPa
            PMPa = P / 1e6;
            T = obj.TempInChamber;
            cCO2 = obj.co2(PMPa);
            cH2O = obj.h2o(PMPa);
            XH2O = obj.Xf(PMPa);
            XH2Od= cH2O ./ (cH2O + cCO2);

            SiO2 = obj.oxides_composition(1);
            [beta, TL, TS] = obj.betaeq(P, XH2Od, T, SiO2);
            rho_m = obj.Rhom(cH2O, cCO2);    % unit density of melt phase
            rho_g = obj.Rhog(PMPa, XH2O, T);    % unit density of gas
            
            t1 = cCO2 + cH2O - 1;
            t9 = rho_m .* (beta - 1) .* (obj.QCO2 .* t1 + obj.QH2O .* t1 + obj.Qc .* (cCO2 + cH2O));
            t12 = beta .* obj.rhoc0 .* (obj.QCO2 + obj.QH2O);
            alpha = max(0,1 ./ (obj.Qc .* rho_g + t12 + t9) .* (t9 + t12));

            rhocFrac = (1.0 - alpha) .* beta .* obj.rhoc0;
            rhogFrac = alpha .* rho_g;
            rhomFrac = (1.0 - alpha) .* (1.0 - beta) .* rho_m;
            rhoFrac = rhogFrac + rhomFrac + rhocFrac;

            obj.rhoInChamber = rhoFrac;
            obj.rhoCrystalInChamber = rhocFrac;
            
            McFrac = rhocFrac / rhoFrac;
            MgFrac = rhogFrac / rhoFrac;
            MmFrac = rhomFrac / rhoFrac;
            obj.uInChamber = McFrac * obj.Cv_c * T +...
                             MmFrac * obj.Cv_m * T +...
                             MgFrac * obj.Ugas(PMPa, XH2O, T);
        end
    end


    methods(Static)
        function [beq, TL, TS] = betaeq(P, XH2Od, T, SiO2)
            P = P / 1e8; % Pa to kbar
            aS=854.0896; bS=6; cS=224.0896; dS=80; eS=0.357; fS=6;
            aL=1205.7; bL=6; cL=285.7; dL=200; eL=0.7; fL=11;
            A=1287.6; B=-20.154;
            
            dTliq=A+B*SiO2;
            TS = 273.15 + aS+bS.*P-XH2Od.*(cS+fS.*P-dS./(P+eS));
            TL = 273.15 + aL+bL*P-XH2Od.*(cL+fL.*P-dL./(P+eL))+dTliq;
            x=(T-TS)./(TL-TS);

            aF=-4.974; bF=28.623; cF=-52.708; dF=34.816;
            beq = (1+exp(aF+bF.*(x)+cF.*(x).^2+dF.*(x).^3)).^(-1);
        end

        function v = Giordano_2008_visc(H2Ot, T, Composition)
            %Get the Giordano2008 VFT parameters from their script
            %Giordano2008_Model.mat
            outvalues = StateSolver.Giordano2008_Model(H2Ot, Composition);
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
            [ncomps, xmf_t] = StateSolver.molepct_grd(wtm);
            
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
            
        function [nox, xmf] = molepct_grd(wtm)
            [nr, nc] = size(wtm);
            % [n x]=MOLEPCT(X) Calculates mole percent oxides from wt %
            % 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1
            % Output: mole fractions of equivalent
            mw = [60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];
            
            %Replicate the molecular weight row by the number of spatial nodes used
            %(number of data rows)
            mw_rep = repmat(mw, nr, 1);
            
            %Perform the calculation using matrix algebra
            wtn = [wtm(:,1:10).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12)) wtm(:,11) 0.5.*wtm(:,12).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12))];
            mp=wtn./mw_rep;
            mpv= 100*(mp./sum(mp,2));
            xmf=mpv;
            
            [~, nc]=size(xmf);
            nox = nc;
        end
            
        function Theta=theta(bx)
            be0=0.45;
            % teta1: experiment fit (different crystallization temperature
            % of oxides in magma)
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


        function etar = bubble_viscosity(alpha)
            % nr = 1.5;
            % b = (1-nr) / (nr * 0.3);
            % etar = 1.0 / (1.0 + b * alpha);
            etar = 1;
        end
    end
end

