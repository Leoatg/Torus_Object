
%LG Pulse

classdef (HandleCompatible) LG_pulse < LG_mode
    
    properties
        force_pulse_time=700;
        pulse_totaltime=70;
        
        ada_coe=0.2;%pulse adiabatic parameter
        ada_shift=30;
        
        pulse_rotate_t=10;%after this time the pulse rotates
    end
    
    methods
        function obj=LG_pulse()
            obj.PL=0.65;
        end
    end
end