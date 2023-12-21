classdef IlsaSettings < matlab.mixin.Copyable
    %ILSASETTINGS Define ILSA settings options
    
    properties
        MASS_TOLERANCE
        LUBRICATION_TOLERANCE
        FRONT_TOLERANCE
        MAX_NONLINEAR_ITERATION
        MAX_LUBRICATION_ITERATION
        MAX_TIP_ELEMENTS
        MIN_WIDTH
        FRONT_RELAXATION
        ASYMPTOTIC_TOLERANCE
        ASYMPTOTIC_MAX_ITERATIONS
        ZERO_FRACTURE_ELEMENTS           % number of elements activated for initial fracture
        FIXED_BOTTOM_TIP                 % if true bottom tip element is fixed
 
        IS_CONSTANT_DENSITY              % if true density is constant
        IS_CONSTANT_TEMPERATURE          % if true energy conservation law is not considered 
        IS_MAGMA_WEIGHT_ON
        IS_CRYSTALLIZATION_ON            % crystallization into energy conservation law
        CRYSTALLIZATION_TYPE
        IS_VISCOUS_DISSIPATION_ON        % viscous heating in energy conseravtion law
        IS_PRESSURE_FORCES_WORK_ON       % accounting pressure forces work in energy conservation law
        
        lubrication_relaxation_iter = [10 20 30 50 200]
        lubrication_relaxation_coef = [1.0 0.7 0.5 0.3 0.1]
        
        IS_HOST_TEMPERATURE_ON
        IS_VARIABLE_HOST_CONDUCTIVITY
        FRACTURE_TYPE = "KGD" % "KGD", "PKN"
        simlog
    end %properties

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = IlsaSettings(settingsInput)
            if nargin > 0
                obj.MASS_TOLERANCE = 1e-4;
                obj.LUBRICATION_TOLERANCE = 1e-6;
                obj.FRONT_TOLERANCE = 1e-3;
                obj.MAX_NONLINEAR_ITERATION = 30;
                obj.MAX_LUBRICATION_ITERATION = 150;
                obj.MAX_TIP_ELEMENTS = 2;
                obj.MIN_WIDTH = settingsInput.MIN_WIDTH;
                obj.FRONT_RELAXATION = 1;
                obj.ASYMPTOTIC_TOLERANCE = 1e-6;
                obj.ASYMPTOTIC_MAX_ITERATIONS = 200;
                obj.ZERO_FRACTURE_ELEMENTS = 5;
                obj.FIXED_BOTTOM_TIP = settingsInput.FIXED_BOTTOM_TIP;

                obj.IS_CONSTANT_DENSITY        = settingsInput.IS_CONSTANT_DENSITY;
                obj.IS_CONSTANT_TEMPERATURE    = settingsInput.IS_CONSTANT_TEMPERATURE;
                obj.IS_MAGMA_WEIGHT_ON         = settingsInput.IS_MAGMA_WEIGHT_ON;
                obj.CRYSTALLIZATION_TYPE       = settingsInput.CRYSTALLIZATION_TYPE;
                obj.IS_CRYSTALLIZATION_ON      = settingsInput.IS_CRYSTALLIZATION_ON;
                obj.IS_VISCOUS_DISSIPATION_ON  = settingsInput.IS_VISCOUS_DISSIPATION_ON;
                obj.IS_PRESSURE_FORCES_WORK_ON = settingsInput.IS_PRESSURE_FORCES_WORK_ON;
                obj.IS_HOST_TEMPERATURE_ON     = settingsInput.IS_HOST_TEMPERATURE_ON;
                obj.IS_VARIABLE_HOST_CONDUCTIVITY = settingsInput.IS_VARIABLE_HOST_CONDUCTIVITY;
            end
        end %default constructor
        
        
        function relaxation_coef = get_lubrication_relaxation_coef(obj, nIter)
            for i = 1:length(obj.lubrication_relaxation_iter)
                if nIter < obj.lubrication_relaxation_iter(i)
                    relaxation_coef = obj.lubrication_relaxation_coef(i);
                    break;
                end
            end
        end %get_lubrication_relaxation_coef
        
        function crystallization_coef = get_crystallization_coef(obj)
            if obj.IS_CRYSTALLIZATION_ON
                crystallization_coef = 1.0;
            else
                crystallization_coef = 0.0;
            end
        end
        
        function pressure_work_coef = get_pressure_work_coef(obj)
            if obj.IS_PRESSURE_FORCES_WORK_ON
                pressure_work_coef = 1.0;
            else
                pressure_work_coef = 0.0;
            end
        end

        function viscous_dissipation_coef = get_viscous_dissipation_coef(obj)
            if obj.IS_VISCOUS_DISSIPATION_ON
                viscous_dissipation_coef = 1.0;
            else
                viscous_dissipation_coef = 0.0;
            end
        end
    end %methods
end