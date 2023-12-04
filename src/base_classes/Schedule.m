classdef Schedule < matlab.mixin.Copyable
    %SCHEDULE Summary of this class goes here
    
    properties
        startTime
        endTime
        timeList
        dtList
        pumpingSchedule
        xPerforation
        rhoPerforation
        dM
    end

    methods(Access = protected)
        function cp = copyElement(obj)
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(obj);
        end
    end
    
    methods
        function obj = Schedule(timeList, dtList, pumpingSchedule, xPerforation, rhoPerforation)
            if nargin > 0
                obj.startTime = timeList(1);
                obj.endTime = timeList(end);
                obj.timeList = timeList;
                obj.dtList = dtList;
                obj.pumpingSchedule = pumpingSchedule;
                obj.xPerforation = xPerforation;
                obj.rhoPerforation = rhoPerforation;
            end
        end %constructor
        
        function Q = injected_volume(obj, startTime, endTime)
            dtStep = endTime - startTime;
            T = obj.pumpingSchedule(1, :);
            T(end+1) = obj.endTime;
            flowrate = obj.pumpingSchedule(2, find(T>startTime, 1) - 1);
            Q = flowrate * dtStep;
        end %injected_volume

        function dM = injected_mass(obj, startTime, endTime)
            dt = endTime - startTime;
            T = obj.pumpingSchedule(1, :);
            T(end+1) = obj.endTime;
            flowrate = obj.pumpingSchedule(2, find(T>startTime, 1) - 1);
            dM = flowrate * dt * obj.rhoPerforation;
        end %injected_mass
        
        function flowrate = get_flowrate(obj, startTime, endTime)
            T = obj.pumpingSchedule(1, :);
            T(end+1) = obj.endTime;
            flowrate = obj.pumpingSchedule(2, find(T>startTime, 1) - 1);
        end %get_flowrate

        function dt = get_dt(obj, startTime)
            dt = obj.dtList(find(obj.timeList > startTime + 0.1, 1) - 1);
        end %get_dt
    end
end

