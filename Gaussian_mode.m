
%LG modes

classdef Gaussian_mode < Grid_basic & GP_parameter & handle
    
    properties
        component_number=2;
        
        Pbar0=3;
        PL=0.4;
        
        r0=0;
        x0=0; %center coordinate
        y0=0;
        Sigma_x=5;
        Sigma_y=5;
        Omega=0;%angular speed of x0 and y0
        u;
        C_phase=0;%constant phase
        
        lambda=0.2815;%760nm
        wave_number_k;
        time=0;
        
        font_size=18;
    end
    
    methods
        
        function obj=Gaussian_mode()
            obj.wave_number_k=2*pi*obj.lambda;
        end
                
        function Build_u(Gobj) 
            if Gobj.r0~=0 && (Gobj.x0~=0 || Gobj.y0~=0)
                if abs(Gobj.r0-sqrt(Gobj.x0^2+Gobj.y0^2))>0.001
                    Gobj.r0=sqrt(Gobj.x0^2+Gobj.y0^2);
                    fprintf('\nWarning: r0 error.\n\n');
                end
            elseif Gobj.r0==0
                if Gobj.x0==0 && Gobj.y0~=0
                    Gobj.r0=Gobj.y0;
                elseif Gobj.x0~=0 && Gobj.y0==0
                    Gobj.r0=Gobj.x0;
                else
                    Gobj.r0=sqrt(Gobj.x0^2+Gobj.y0^2);
                end
            end
            if isempty(Gobj.u)
                [Gobj.x,Gobj.y] = meshgrid(-Gobj.XYmax:Gobj.hspace:Gobj.XYmax);
                Gobj.u=sqrt(Gobj.P0*Gobj.Pbar0)*exp(-((Gobj.x-Gobj.x0).^2/(2*Gobj.Sigma_x^2)+(Gobj.y-Gobj.y0).^2/(2*Gobj.Sigma_y^2)));
                Gobj.u=Gobj.u*exp(1i*Gobj.C_phase);%add a constant phase
            elseif isempty(Gobj.x)
                [Gobj.x,Gobj.y] = meshgrid(-Gobj.XYmax:Gobj.hspace:Gobj.XYmax);
            else
                fprintf('\nWarning: u has been defined.\n\n');
            end
        end
        
        function Equal=Check_Parameter(LGobj1,LGobj2)% 2 is the 'old' one
            if LGobj2.component_number==1
                Equal= (LGobj1.Pbar0==LGobj2.Pbar0) ...
                    && (LGobj1.Sigma_x==LGobj2.Sigma_x) && (LGobj1.Sigma_y==LGobj2.Sigma_y)...
                    && (LGobj1.x0==LGobj2.x0) && (LGobj1.y0==LGobj2.y0)...
                    && (LGobj1.ua==LGobj2.ua) ...
                    && (LGobj1.R==LGobj2.R) && (LGobj1.J==LGobj2.J)...
                    && (LGobj1.gR==LGobj2.gR) && (LGobj1.GammaR==LGobj2.GammaR)...
                    && (LGobj1.component_number==LGobj2.component_number)...
                    && (LGobj1.XYmax==LGobj2.XYmax) && (LGobj1.dt==LGobj2.dt);
            elseif LGobj2.component_number==2
                Equal= (LGobj1.Pbar0==LGobj2.Pbar0) && (LGobj1.PL==LGobj2.PL)  ...
                    && (LGobj1.Sigma_x==LGobj2.Sigma_x) && (LGobj1.Sigma_y==LGobj2.Sigma_y)...
                    && (LGobj1.x0==LGobj2.x0) && (LGobj1.y0==LGobj2.y0)...
                    && (LGobj1.ua==LGobj2.ua) ...
                    && (LGobj1.R==LGobj2.R) && (LGobj1.J==LGobj2.J)...
                    && (LGobj1.gR==LGobj2.gR) && (LGobj1.GammaR==LGobj2.GammaR)...
                    && (LGobj1.component_number==LGobj2.component_number)...
                    && (LGobj1.XYmax==LGobj2.XYmax) && (LGobj1.dt==LGobj2.dt);
            else
                fprintf('\nError: component number incorrect.\n\n');
                return;
            end
        end %function Check_Parameter ends
        
       
        function Update_u(Gobj)
            if isempty(Gobj.u) || isempty(Gobj.x)
                [Gobj.x,Gobj.y] = meshgrid(-Gobj.XYmax:Gobj.hspace:Gobj.XYmax);
                Gobj.u=sqrt(Gobj.P0*Gobj.Pbar0)*exp(-((Gobj.x-Gobj.x0).^2/(2*Gobj.Sigma_x^2)+(Gobj.y-Gobj.y0).^2/(2*Gobj.Sigma_y^2)));
                fprintf('\nWarning: u has not been defined.\n\n');
            else
                if Gobj.r0~=0
                    Gobj.x0=cos(double(Gobj.Omega*Gobj.time))*Gobj.r0;
                    Gobj.y0=sin(double(Gobj.Omega*Gobj.time))*Gobj.r0;
                    Gobj.u=sqrt(Gobj.P0*Gobj.Pbar0)*exp(-((Gobj.x-Gobj.x0).^2/(2*Gobj.Sigma_x^2)+(Gobj.y-Gobj.y0).^2/(2*Gobj.Sigma_y^2)));
                elseif Gobj.r0==0
                    %do nothing
                end
            end
        end
        
    end
    
    %mesh
    methods
        function [handle_fig,handle_mesh]=mesh_u(Gobj)
               if ~isempty(Gobj.u)
                  handle_fig=figure('position',[500 400 700 700],'renderer','opengl');
                  handle_mesh=mesh(Gobj.x,Gobj.y,abs(Gobj.u).^2);
                  set(gca,'fontsize',Gobj.font_size);
                  axis square;
                  colorbar;
                  view(0,90);
                  xlabel('x');
                  ylabel('y');
               else
                   fprintf('\nError: u is not yet built.\n\n');
               end 
        end
           
        
        function clear_xy(Gobj)
            Gobj.x=[];
            Gobj.y=[];
        end
        
        function clear_all(Gobj)
            Gobj.x=[];
            Gobj.y=[];
            Gobj.u=[];
            Gobj.time=0;
        end
    end
end












