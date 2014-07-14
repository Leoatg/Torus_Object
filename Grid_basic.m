

%Grid basic

classdef Grid_basic < handle
    
    properties
        N=1024;
        XYmax=60;
        hspace=0.1173;
        dt;
        TotalTime;
        x;
        y;
        damping_shift; %damping distance from the boundary
        a_damping; %slope for the damping mesa
    end
    
    methods
        function obj=Grid_basic()
            obj.dt=floor((obj.hspace^2/pi)*1e3)*1e-3;
        end
        
        function copy_Grid(obj,Gridobj)% right to left
            obj.N=Gridobj.N;
            obj.XYmax=Gridobj.XYmax;
            obj.hspace=Gridobj.hspace;
            obj.dt=Gridobj.dt;
            obj.damping_shift=Gridobj.damping_shift;
            obj.a_damping=Gridobj.a_damping;
            if ~isempty(Gridobj.TotalTime)
                obj.TotalTime=Gridobj.TotalTime;
            else
                fprintf('\nWarning: TotalTime is not defined\n\n');
            end
        end
    end
    
    
end