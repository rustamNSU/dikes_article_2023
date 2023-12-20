classdef Reservoir < matlab.mixin.Copyable
    %Reservoir contain reservoir properties
    properties
        E
        nu
        toughness
        leakoff
        Ep
        Kp
        Cp % Leakoff coefficient (zero), need to asymptotic algorithm for tip elements
        elasticityMatrix
        sigmaH
        temperatureH
        gravity = 9.81
        alpha
        
        k
        cp
        rho
        Ly
        Ny
    end

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = Reservoir(input, settings)
            if nargin > 0    
                obj.E  = input.Young;
                obj.nu = input.nu;
                obj.toughness = input.KIc;
                obj.rho = input.rho;
                obj.gravity = input.gravity;
                obj.sigmaH = input.sigmaH;

                obj.Ep = obj.E / (1.0 - obj.nu^2);
                obj.Kp = sqrt(32.0 / pi) * obj.toughness;
                obj.leakoff = 1e-20;
                obj.Cp = 1e-20;
                
                if (settings.IS_CONSTANT_TEMPERATURE == false)
                    obj.cp = input.cp;
                    obj.k = input.k;
                    obj.temperatureH = input.temperatureH;
                    obj.Ly = input.Ly;
                    obj.Ny = input.Ny;
                end
            end
        end
    end
end

