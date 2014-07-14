

%GP equation parameters

classdef GP_parameter < handle
    
    properties
        
        ua=0.0077; %same spin polariton-polariton scattering
        GammaR_ratio=1.5;
        GammaC=1;  % GammaR_ratio=GammaR/GammaC
        R=0.0084; %stimulated scattering rate
        
        J=0.5;
        
        ub;
        gR;
        GammaR;
        P0; %threshold pumping power Pth
        
        A_Omega=0; %angular momentum in rotating frame
    end
    
    methods
        
        function obj=GP_parameter()
            obj.ub=-0.1*obj.ua; %cross spin scattering
            obj.gR=2*obj.ua; %LP interaction with incoherent reservoir
            obj.GammaR=obj.GammaR_ratio*obj.GammaC;
            obj.P0=(obj.GammaR*obj.GammaC)/obj.R;
        end
        
        function Update_GP_parameter(GPobj)
            GPobj.GammaR=GPobj.GammaR_ratio*GPobj.GammaC;
            GPobj.P0=(GPobj.GammaR*GPobj.GammaC)/GPobj.R;
        end
        
        
        function copy_GP(obj,GPobj) %right to left
            obj.ua=GPobj.ua;
            obj.R=GPobj.R;
            obj.J=GPobj.J;
            obj.ub=GPobj.ub;
            obj.gR=GPobj.gR;
            obj.GammaR_ratio=GPobj.GammaR_ratio;
            obj.GammaC=GPobj.GammaC;
            obj.GammaR=GPobj.GammaR;
            obj.P0=GPobj.P0;
        end
    end
   
end