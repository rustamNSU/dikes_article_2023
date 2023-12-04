classdef FractureElements < matlab.mixin.Copyable
    %FRACTUREELEMENTS Summary of this class goes here
    
    properties
        width          % [m]
        pressure       % [Pa]
        rho            % [kg/m^3]
        temperature    % [K]
        u              % [J/kg]
        mu             % [Pa*s]
        rockTemperature % class RockTemperature
        crystallization % class CrystallizationSolver
        
        alpha
        beta
        betaeq
        rhoc   % crystal density in dike (relative to total volume)
        rhom   % melt phase density in dike (relative to total volume)
        rhog   % gas density in dike (relative to total volume)
        ch2o 
        cco2
        xh2o
        xh2od
        Mc     % crystal mass content
        Mg     % gas mass content
        Mm     % melt mass content
        TL     % crystal liqidus temperature
        TS     % crystal solidus temperature
        
        flux_velocity
        mass_rate
        front_velocity

        channel
        tip
        survey
        distanceToSurvey
        leftTip
        rightTip
        mesh % class Mesh
        time
        perforation
        front
    end %properties

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
            
            % Deep copy of handle properties
            if obj.rockTemperature ~= 0
                cp.rockTemperature = copy(obj.rockTemperature);
            end
            if obj.crystallization ~= 0
                cp.crystallization = copy(obj.crystallization);
            end
        end
    end
    
    methods
        function obj = FractureElements(mesh)
            if nargin > 0
                n = mesh.n;
                obj.mesh = mesh;
                obj.width = zeros(1, n);
                obj.pressure = zeros(1, n);
                obj.rho = zeros(1, n);
                obj.temperature = zeros(1, n);
                obj.u             = zeros(1, n);
                obj.mu = zeros(1, n);
                obj.flux_velocity = zeros(1, n+1);
                obj.mass_rate = zeros(1, n+1);

                obj.alpha  = zeros(1, n);
                obj.beta   = zeros(1, n);
                obj.betaeq = zeros(1, n);
                obj.rhoc   = zeros(1, n);
                obj.rhom   = zeros(1, n);
                obj.rhog   = zeros(1, n);
                obj.ch2o   = zeros(1, n);
                obj.cco2   = zeros(1, n);
                obj.xh2o   = zeros(1, n);
                obj.xh2od  = zeros(1, n);
                obj.Mc     = zeros(1, n);
                obj.Mg     = zeros(1, n);
                obj.Mm     = zeros(1, n);
                obj.TL     = zeros(1, n);
                obj.TS     = zeros(1, n);
            end
        end %default constructor
        
        
        function obj = initialize_zero_fracture(obj, settings, xPerforation, startTime)
            obj.time = startTime;
            hasPerforation = (obj.mesh.leftBoundary < xPerforation &...
                obj.mesh.rightBoundary >= xPerforation);
            obj.perforation = find(hasPerforation);
            % obj.channel = obj.perforation(1) - 1:1:obj.perforation(1) + 1;
            obj.channel = obj.perforation;
            obj.survey = [obj.channel(1), obj.channel(end)];
            obj.tip = obj.survey + [-1, 1];
%             obj.distanceToSurvey = (obj.mesh.dx(obj.survey) + obj.mesh.dx(obj.tip)) / 2.0;
            obj.distanceToSurvey = 0.51 * obj.mesh.dx(obj.survey);
            obj.front = obj.mesh.xc(obj.survey) + [-1, 1].*obj.distanceToSurvey;
            obj.leftTip = obj.tip(1);
            obj.rightTip = obj.tip(end);
            obj.width([obj.channel obj.tip]) = settings.MIN_WIDTH;
            obj.front_velocity = [0.0, 0.0];
        end %initialize_zero_fracture


        function obj = initialize_rock_temperature(obj, Ly, Ny, temperatureH)
            obj.rockTemperature = RockTemperature(Ly, obj.mesh.n, Ny, temperatureH)
        end %initialize_rock_temperature
        
        
        function elements = get_fracture_elements(obj)
            elements = sort([obj.channel, obj.tip]);
        end %get_fracture_elements
        
        
        function V = get_fracture_volume(obj)
            V = dot(obj.width, obj.mesh.dx);
        end %get_fracture_volume
        
        
        function m = get_fracture_fluid_mass(obj)
            mElements = obj.rho .* obj.width;
            m = dot(mElements, obj.mesh.dx);
        end %get_fracture_fluid_mass
        
        
        function m = get_mass(obj)
            mElements = obj.rho .* obj.width;
            m = mElements .* obj.mesh.dx;
        end
        
        
        function width = get_fracture_width(obj)
            width = obj.width(sort([obj.channel, obj.tip]));
        end
        
        
        function covering = get_elements_covering(obj)
            n = obj.mesh.n;
            covering = zeros(1, n);
            boundary = obj.mesh.get_boundary((1:n));
            lmin = obj.front(1);
            lmax = obj.front(2);
            
            for i = 1:n
                lb = boundary(i);
                rb = boundary(i+1);
                dx = obj.mesh.dx(i);
                
                if (lb > lmin && rb < lmax)
                    covering(i) = 1.0;
                elseif (rb < lmin || lb > lmax)
                    covering(i) = 0.0;
                elseif (lb < lmin)
                    covering(i) = (rb - lmin) / dx;
                elseif (rb > lmax)
                    covering(i) = (lmax - lb) / dx;
                else
                    error('Stop calculation because bug in FractureElements.get_elements_covering()');
                end
            end
        end
        
        
        function set_constant_temperature(obj, T0)
            obj.temperature(:) = T0;
        end %set_constant_temperature
        

        function set_channel(obj, wNewChannel, pressureChannel)
            isEqualLength = (length(wNewChannel) == length(obj.channel)) &&...
                (length(pressureChannel) == length(obj.channel));
            if ~isEqualLength
                error('error in set_channel: incorrect input argument');
            end
            obj.width(obj.channel) = wNewChannel;
            obj.pressure(obj.channel) = pressureChannel;
        end %set_channel

        
        function set_tip_pressure(obj, pressureTip)
            isEqualLength = length(pressureTip) == length(obj.tip);
            if ~isEqualLength
                error('error in set_tip_pressure: not equal length of arrays');
            end
            obj.pressure(obj.tip) = pressureTip;
        end %set_tip_pressure
        
        
        function set_survey_distance(obj, sNew)
            if not(length(sNew)==length(obj.distanceToSurvey))
                error('error in set_survey_distance(obj, sNew): not equal length of arrays');
            end
            obj.distanceToSurvey = sNew;
            obj.front = obj.mesh.xc(obj.survey) - [1 -1].*sNew;
            
            %% Update tip elements
            leftWing = 1:obj.survey(1)-1;
            rightWing = obj.survey(2)+1:obj.mesh.n;
            obj.leftTip = leftWing(obj.mesh.rightBoundary(leftWing) > obj.front(1));
            obj.rightTip = rightWing(obj.mesh.leftBoundary(rightWing) < obj.front(2));
            obj.tip = [obj.leftTip obj.rightTip];
            obj.width(obj.tip) = 0;
            obj.pressure(obj.tip) = 0;
        end %set_survey_distance
        
        
        %% Replace survey elements to inside tip
        function update_survey_elements(obj)
            obj.channel = [obj.leftTip obj.channel obj.rightTip];
            obj.channel([1 end]) = [];
            obj.leftTip = obj.leftTip(1);
            obj.rightTip = obj.rightTip(end);
            obj.tip = [obj.leftTip obj.rightTip];
            obj.survey = obj.channel([1 end]);
            obj.distanceToSurvey = abs(obj.mesh.xc(obj.survey) - obj.front);
        end %update_survey_elements


        %% Update tip parameters after adding new tip elements
        function update_tip_parameters(obj)
            leftSurvey = obj.survey(1);
            rightSurvey = obj.survey(end);
            leftTip = obj.leftTip;
            rightTip = obj.rightTip;

            obj.rho(leftTip)  = obj.rho(leftSurvey);
            obj.rho(rightTip) = obj.rho(rightSurvey);

            obj.rhoc(leftTip)  = obj.rhoc(leftSurvey);
            obj.rhoc(rightTip) = obj.rhoc(rightSurvey);

            obj.pressure(leftTip)  = obj.pressure(leftSurvey);
            obj.pressure(rightTip) = obj.pressure(rightSurvey);

            obj.temperature(leftTip)  = obj.temperature(leftSurvey);
            obj.temperature(rightTip) = obj.temperature(rightSurvey);

            obj.alpha(leftTip)  = obj.alpha(leftSurvey);
            obj.alpha(rightTip) = obj.alpha(rightSurvey);

            obj.beta(leftTip)  = obj.beta(leftSurvey);
            obj.beta(rightTip) = obj.beta(rightSurvey);

            obj.betaeq(leftTip)  = obj.betaeq(leftSurvey);
            obj.betaeq(rightTip) = obj.betaeq(rightSurvey);

            obj.mu(leftTip)  = obj.mu(leftSurvey);
            obj.mu(rightTip) = obj.mu(rightSurvey);

            obj.Mc(leftTip)  = obj.Mc(leftSurvey);
            obj.Mc(rightTip) = obj.Mc(rightSurvey);

            obj.Mm(leftTip)  = obj.Mm(leftSurvey);
            obj.Mm(rightTip) = obj.Mm(rightSurvey);

            obj.Mg(leftTip)  = obj.Mg(leftSurvey);
            obj.Mg(rightTip) = obj.Mg(rightSurvey);
            if isprop(obj, 'crystallization') && obj.crystallization ~= 0
                if isprop(obj.crystallization, 'csd')
                    for ix = leftTip
                        obj.crystallization.csd(ix, 1:end)  = obj.crystallization.csd(leftSurvey, 1:end);
                    end
                    for ix = rightTip
                        obj.crystallization.csd(ix, 1:end) = obj.crystallization.csd(rightSurvey, 1:end);
                    end
                else
                    for ix = leftTip
                        for i = 1:length(obj.crystallization.betaeqInterp)
                            obj.crystallization.crystallizationSolvers{i}.csd(ix, 1:end)  = obj.crystallization.crystallizationSolvers{i}.csd(leftSurvey, 1:end);
                        end
                    end
                    for ix = rightTip
                        for i = 1:length(obj.crystallization.betaeqInterp)
                            obj.crystallization.crystallizationSolvers{i}.csd(ix, 1:end)  = obj.crystallization.crystallizationSolvers{i}.csd(rightSurvey, 1:end);
                        end
                    end
                end
            end
        end %update_tip_parameters


        function deepCopy(obj, orgObj)
            pl = properties(orgObj);
            for k = 1:length(pl)
               if isprop(obj,pl{k})
                  obj.(pl{k}) = orgObj.(pl{k});
               end
            end
            if isprop(obj, 'rockTemperature') && obj.rockTemperature ~= 0
                obj.rockTemperature = copy(orgObj.rockTemperature);
            end
            if isprop(obj, 'crystallization') && obj.crystallization ~= 0
                obj.crystallization = copy(orgObj.crystallization);
            end
        end
    end %methods
end

