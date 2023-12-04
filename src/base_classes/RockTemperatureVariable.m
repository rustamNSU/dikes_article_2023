classdef RockTemperatureVariable < matlab.mixin.Copyable
    properties
        Ly
        ny
        nx
        yc
        yBoundary
        dy

        channel
        fracElements
        Tm
        
        k0
        k
        rho
        C

        T
        TInf
        penalty = 1;
    end %properties

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = RockTemperatureVariable(reservoir)
            if nargin > 0
                ny = reservoir.Ny;
                Ly = reservoir.Ly;
                nx = length(reservoir.temperatureH);

                obj.yc = zeros(1, ny);
                obj.yBoundary = zeros(1, ny);
                obj.dy = zeros(1, ny);
                obj.Ly = Ly;
                obj.ny = ny;
                obj.nx = nx;
                obj.rho = reservoir.rho;
                obj.C = reservoir.cp;
                obj.k0 = reservoir.k;
                obj.T = zeros(nx, ny);
                obj.TInf = reservoir.temperatureH;
                obj.Tm = reservoir.temperatureH;
                for i = 1:nx
                    obj.T(i, 1:end) = obj.TInf(i);
                end
                obj.create_mesh_pow(2);
                obj.k = obj.updateConductivity(obj.T, obj.k0);
            end
        end %default constructor
        

        function create_mesh(obj)
            ny = obj.ny;
            Ly = obj.Ly;
            yc = zeros(1, ny);
            yb = zeros(1, ny);
            dy = zeros(1, ny);
            for i = 1:ny
                yb(i) = Ly / ny * i;
            end

            yc(1) = 0.5 * yb(1);
            dy(1) = yb(1);
            for i = 2:ny
                yc(i) = 0.5 * (yb(i-1) + yb(i));
                dy(i) = yb(i) - yb(i-1);
            end

            obj.yc = yc;
            obj.yBoundary = yb;
            obj.dy = dy;
        end %create_mesh

        function create_mesh_pow(obj, pow)
            ny = obj.ny;
            Ly = obj.Ly;
            yc = zeros(1, ny);
            yb = zeros(1, ny);
            dy = zeros(1, ny);
            for i = 1:ny
                yb(i) = Ly * (i / ny)^pow;
            end

            yc(1) = 0.5 * yb(1);
            dy(1) = yb(1);
            for i = 2:ny
                yc(i) = 0.5 * (yb(i-1) + yb(i));
                dy(i) = yb(i) - yb(i-1);
            end

            obj.yc = yc;
            obj.yBoundary = yb;
            obj.dy = dy;
        end %create_mesh


        function qY = get_heat_flow(obj)
            dT = obj.T(:, 1) - obj.Tm';
            dy = obj.dy(1) / 2;
            qY = -obj.k(:, 1) .* dT / dy;
        end %get_heat_flow


        function solve_heat_equation(obj, oldObj, dt)
            TInf = obj.TInf;
            dy   = obj.dy;
            yc   = obj.yc;
            ny   = obj.ny;
            rho  = obj.rho;
            C    = obj.C;
            obj.k = obj.updateConductivity(obj.T, obj.k0);
            k = obj.k;
        
            n = 3*ny - 2;
            I = zeros(1, n);
            J = zeros(1, n);
            V = zeros(1, n);
            b = zeros(1, ny);
            for ix = obj.channel
                aR  = 0.5 * (k(ix, 2) + k(ix, 1)) / dy(1) / (yc(2)-yc(1));
                aLB = k(ix, 1) / dy(1) / dy(1) / 2;
                I(1) = 1;
                J(1) = 1;
                V(1) = rho*C/dt + aR + aLB;

                I(2) = 1;
                J(2) = 2;
                V(2) = -aR;
                b(1) = obj.Tm(ix)*aLB + rho*C/dt*oldObj.T(ix, 1);
                ind = 3;
                for iy = 2:ny-1
                    aL = 0.5 * (k(ix, iy) + k(ix, iy-1)) / dy(iy) / (yc(iy)-yc(iy-1));
                    aR = 0.5 * (k(ix, iy+1) + k(ix, iy)) / dy(iy) / (yc(iy+1)-yc(iy));
                    I(ind) = iy;
                    J(ind) = iy-1;
                    V(ind) = -aL;
                    ind = ind + 1;
        
                    I(ind) = iy;
                    J(ind) = iy;
                    V(ind) = rho*C/dt + aL + aR;
                    ind = ind + 1;
        
                    I(ind) = iy;
                    J(ind) = iy+1;
                    V(ind) = -aR;
                    ind = ind + 1;
        
                    b(iy) = rho*C/dt*oldObj.T(ix, iy);
                end
                aL  = 0.5 * (k(ix, ny) + k(ix, ny-1)) / dy(ny) / (yc(ny)-yc(ny-1));
                aRB = k(ix, ny) / dy(ny) / dy(ny) / 2;
                I(ind) = ny;
                J(ind) = ny-1;
                V(ind) = -aL;
                ind = ind + 1;

                I(ind) = ny;
                J(ind) = ny;
                V(ind) = rho*C/dt + aL + aRB;
                b(ny) = obj.TInf(ix)*aRB + rho*C/dt*oldObj.T(ix, ny);
        
                A = sparse(I, J, V, ny, ny);
                x = A\b';
                obj.T(ix,:) = x';
            end
        end
    end %methods


    methods(Static)
        function k = updateConductivity(T, k0)
            b = 2e-3;
            k = k0 ./ (1.0 + b * (T-273.15));
        end
    end
end

