
%LG Pulse

classdef (HandleCompatible) Gaussian_pulse < Gaussian_mode
    
    properties
        force_pulse_time=700;
        pulse_totaltime=70;
        
        ada_coe=0.2;%pulse adiabatic parameter
        ada_shift=30;
    end
    
    methods
        function obj=Gaussian_pulse()
            obj.PL=0.65;
        end
    end
end