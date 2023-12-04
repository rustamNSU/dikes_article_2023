classdef Mesh < matlab.mixin.Copyable
    %MESH 
    
    properties
        xc
        leftBoundary
        rightBoundary
        dx
        n
    end %properties

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = Mesh(xc, dx)
            if nargin > 0
                if not (length(xc) == length(dx))
                    error('Mesh created error, different length of xc and dx')
                end
                obj.xc = xc;
                obj.dx = dx;
                obj.n = length(xc);
                obj.leftBoundary = xc - dx./2;
                obj.rightBoundary = xc + dx./2;
            end
        end

        function boundary = get_boundary(obj, elements)
            boundary = obj.leftBoundary(elements);
            boundary = [boundary, obj.rightBoundary(elements(end))];
        end
    end %methods
end

