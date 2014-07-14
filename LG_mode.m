
%LG modes

classdef LG_mode < Grid_basic & GP_parameter & handle
    
    properties
        component_number=2;
        
        Pbar0=3;
        PL=0.4;
        
        l=5;
        p=8;
        w0=24;
        zR=5;
        z=0;
        Omega=0;
        lambda=0.2815;%760nm
        wave_number_k;
        time=0;
        
        %build u parameters
        phi;r;u;
        u0;%initial u without phase evolution
        
        font_size=18;
    end
    
    methods
        
        function obj=LG_mode()
            obj.wave_number_k=2*pi*obj.lambda;
        end
        
        function Build_u(LGobj)
            if isempty(LGobj.u)
                [LGobj.x,LGobj.y] = meshgrid(-LGobj.XYmax:LGobj.hspace:LGobj.XYmax);
                LGobj.r=sqrt(LGobj.x.^2+LGobj.y.^2);
                %LGobj.phi=atan2(LGobj.y,LGobj.x);
                if LGobj.z==0
                    LGobj.u=double(exp(1i*LGobj.Omega*LGobj.time))*(1/LGobj.w(LGobj.z))*((LGobj.r*sqrt(2)/LGobj.w(LGobj.z)).^abs(LGobj.l)).* exp(-LGobj.r.^2/LGobj.w(LGobj.z)^2).* abs(polyval(LaguerreGen(LGobj.p,abs(LGobj.l)), (2*LGobj.r.^2)/LGobj.w(LGobj.z)^2))...
                        .*exp(1i*LGobj.l*atan2(LGobj.y,LGobj.x)).*exp(1i*(2*LGobj.p+abs(LGobj.l)+1)*LGobj.Zeta(LGobj.z));
                    
                else
                    LGobj.u=double(exp(1i*LGobj.Omega*LGobj.time))*(1/LGobj.w(LGobj.z))*((LGobj.r*sqrt(2)/LGobj.w(LGobj.z)).^abs(LGobj.l)).* exp(-LGobj.r.^2/LGobj.w(LGobj.z)^2).*abs(polyval(LaguerreGen(LGobj.p,abs(LGobj.l)), (2*LGobj.r.^2)/LGobj.w(LGobj.z)^2))...
                        .*exp(1i*LGobj.l*atan2(LGobj.y,LGobj.x)).*exp(1i*(2*LGobj.p+abs(LGobj.l)+1)*LGobj.Zeta(LGobj.z)).*exp(1i*LGobj.wave_number_k*LGobj.r.^2/(2*LGobj.R_z(LGobj.z)));
                end
                %normalized by the square
                P_max=max(max(abs(LGobj.u)));
                LGobj.u=sqrt(LGobj.P0*LGobj.Pbar0)*LGobj.u/P_max;
                LGobj.u0=LGobj.u;
            else
                fprintf('\nError: u has been defined.\n\n');
            end
        end
        
        function z_o=w(LGobj,z_i)
            z_o=LGobj.w0*sqrt(1+(z_i/LGobj.zR)^2);
        end
        
        function z_o=R_z(LGobj,z_i)
            if z_i~=0
                z_o=z_i*(1+(LGobj.zR/z_i)^2);
            else
                fprintf('Error: z = 0.\n')
                return;
            end
        end
        
        function z_o=Zeta(LGobj,z_i)
            z_o=atan(z_i/LGobj.zR);
        end
        
        function Update_u(LGobj)
            if isempty(LGobj.u)
                LGobj.Build_u();
                fprintf('\nWarning: LG mode not exist.\n\n');
            end
            LGobj.u=exp(1i*LGobj.Omega*LGobj.time)*LGobj.u0;
        end
        
    end
    
    %mesh
    methods
        function [handle_fig,handle_mesh]=mesh_u(LGobj)
            if ~isempty(LGobj.u)
                handle_fig=figure('position',[500 400 700 700],'renderer','opengl');
                handle_mesh=mesh(LGobj.x,LGobj.y,abs(LGobj.u).^2);
                set(gca,'fontsize',LGobj.font_size);
                axis square;
                colorbar;
                view(0,90);
                xlabel('x');
                ylabel('y');
            else
                fprintf('\nError: u is not yet built.\n\n');
            end
        end
        
        function Equal=Check_Parameter(LGobj1,LGobj2)% 2 is the 'old' one
            if LGobj2.component_number==1
                Equal= (LGobj1.Pbar0==LGobj2.Pbar0) && (LGobj1.dt==LGobj2.dt)...
                    && (LGobj1.XYmax==LGobj2.XYmax) && (LGobj1.N==LGobj2.N)...
                    && (LGobj1.l==LGobj2.l) && (LGobj1.p==LGobj2.p)...
                    && (LGobj1.w0==LGobj2.w0) && (LGobj1.z==LGobj2.z)...
                    && (LGobj1.ua==LGobj2.ua) ...
                    && (LGobj1.R==LGobj2.R) && (LGobj1.J==LGobj2.J)...
                    && (LGobj1.gR==LGobj2.gR) && (LGobj1.GammaR==LGobj2.GammaR)...
                    && (LGobj1.component_number==LGobj2.component_number)...
                    && (LGobj1.damping_shift==LGobj2.damping_shift)...
                    && (LGobj1.a_damping==LGobj2.a_damping);
            elseif LGobj1.component_number==2
                Equal= (LGobj1.Pbar0==LGobj2.Pbar0) && (LGobj1.PL==LGobj2.PL)...
                    && (LGobj1.dt==LGobj2.dt)...
                    && (LGobj1.XYmax==LGobj2.XYmax) && (LGobj1.N==LGobj2.N)...
                    && (LGobj1.l==LGobj2.l) && (LGobj1.p==LGobj2.p)...
                    && (LGobj1.w0==LGobj2.w0) && (LGobj1.z==LGobj2.z)...
                    && (LGobj1.ua==LGobj2.ua) && (LGobj1.ub==LGobj2.ub)...
                    && (LGobj1.R==LGobj2.R) && (LGobj1.J==LGobj2.J)...
                    && (LGobj1.gR==LGobj2.gR) && (LGobj1.GammaR==LGobj2.GammaR)...
                    && (LGobj1.component_number==LGobj2.component_number);
            else
                fprintf('\nError: component number incorrect.\n\n');
                return;
            end
        end %function Check_Parameter ends
        
        function clear_xy(LGobj)
            LGobj.x=[];
            LGobj.y=[];
        end
        
        function clear_all(LGobj)
            LGobj.x=[];
            LGobj.y=[];
            LGobj.r=[];
            LGobj.phi=[];
            LGobj.u=[];
            LGobj.u0=[];
            LGobj.time=0;
        end
    end
end












