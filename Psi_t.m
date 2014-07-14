

%psi class


classdef (HandleCompatible) Psi_t < GP_parameter & Grid_basic
    
    properties
        component_number=2;
        
        psi;phi;%single component case
        n;
        psi_trim;%record psi
        n_trim;%record reservoir
       
        psiL; psiR;%two component case
        nL;   nR;
        phiL; phiR;
        psiL_trim;%record parameter
        psiR_trim;
        nL_trim;
        nR_trim;
        
        Dk=0;%TE-TM splitting
        
        nc; %Strokes parameter
        sx;sy;sz;
        theta;
        record_sxyz_int=0;
        sx_int;sy_int;sz_int;
        
        %initial state parameter
        Ini_m=0;%angular momemtum for initial state
        Ini_kind=2; %1: Gaussian, 2: Torodial
        IR=0.1;
        IR_zmax=0.2;
        IR_sigma_x=10; %parameters for Gaussian
        IR_sigma_y=10;
        Ini_r=4;  %parameters for torodial %cross section width
        Ini_R=30; %bulk radius
        
        
        %record trimmed psi parameters
        record_psi=0;
        re_psi_begin=70;
        re_psi_end;
        re_psi_dt=1;
        
        step_record_psi=1;
        record_psi_t;
        record_finial_state=0;
        
        %record Time and Energy parameters
        Time;
        Energy;
        Energy_unmal;
        Lz;
        Psi_int;
        mod_psiLR_square;
        E_t_interval=0.1
        E_step_interval;
        step_time=1;
        step_E=1;
        step_i;
         
        %integration area
        fixed_int_r=1; %0: variable, 1:fixed
        int_tolerant=1e-3;
        int_r=40;
        int_area;
        
        TotalStep;
        font_size=18;
        zoom_factor=2;
    end
    
    methods
       
        function initialize(psiobj)
            psiobj.TotalStep=floor(psiobj.TotalTime/psiobj.dt);
            if psiobj.component_number==1 %single component case
                %initialize psi and n--Single
                if isempty(psiobj.psi)
                    if isempty(psiobj.x)
                        [psiobj.x,psiobj.y] = meshgrid(-psiobj.XYmax:psiobj.hspace:psiobj.XYmax);
                    end
                    if psiobj.Ini_kind==1 %Gaussian initial state
                    psiobj.psi=psiobj.IR_zmax*exp(-(psiobj.x.^2/psiobj.IR_sigma_x^2+psiobj.y.^2/psiobj.IR_sigma_y^2));
                    elseif psiobj.Ini_kind==2
                    psiobj.psi=psiobj.IR_zmax*exp(-(sqrt(psiobj.x.^2+psiobj.y.^2)-psiobj.Ini_R).^2/(2*psiobj.Ini_r^2))/(psiobj.Ini_r*sqrt(2*pi));
                    end
                    psiobj.psi=psiobj.psi.*double(exp(1i*psiobj.Ini_m*atan2(psiobj.y,psiobj.x)));%initial angular momentum
                    psiobj.n=abs(psiobj.psi);
                elseif isempty(psiobj.x)
                    [psiobj.x,psiobj.y] = meshgrid(-psiobj.XYmax:psiobj.hspace:psiobj.XYmax);
                else
                    fprintf('\nWarning: psi has been initialized.\n\n');
                end
                
                %initialize record psi time
                if psiobj.record_psi~=0 
                    if  isempty(psiobj.TotalTime)
                        psiobj.re_psi_end=psiobj.TotalTime;
                        fprintf('\nWarning: TotalTime is not defined.\n\n');
                    else
                        psiobj.record_psi_t=psiobj.re_psi_begin:psiobj.re_psi_dt:psiobj.re_psi_end;
                        psiobj.psi_trim=complex(zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t)),0);
                        psiobj.n_trim=zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t));
                        psiobj.TotalStep=floor(psiobj.TotalTime/psiobj.dt);
                    end
                end
                
            elseif psiobj.component_number==2 %two component case
                               
                %initialize psiLR and nLR
                if isempty(psiobj.psiL)
                    if isempty(psiobj.x)
                        [psiobj.x,psiobj.y] = meshgrid(-psiobj.XYmax:psiobj.hspace:psiobj.XYmax);
                    end
                    if psiobj.Ini_kind==1 %Gaussian initial state
                        psiobj.psiR=psiobj.IR_zmax*exp(-(psiobj.x.^2/psiobj.IR_sigma_x^2+psiobj.y.^2/psiobj.IR_sigma_y^2));
                    elseif psiobj.Ini_kind==2
                        psiobj.psiR=psiobj.IR_zmax*exp(-(sqrt(psiobj.x.^2+psiobj.y.^2)-psiobj.Ini_R).^2/(2*psiobj.Ini_r^2))/(psiobj.Ini_r*sqrt(2*pi));
                    end
                    psiobj.psiR=psiobj.psiR.*double(exp(1i*psiobj.Ini_m*atan2(psiobj.y,psiobj.x)));%initial angular momentum
                    psiobj.psiL=psiobj.IR*psiobj.psiR;
                    psiobj.nL=abs(psiobj.psiL);
                    psiobj.nR=abs(psiobj.psiR);
                elseif isempty(psiobj.x)
                    [psiobj.x,psiobj.y] = meshgrid(-psiobj.XYmax:psiobj.hspace:psiobj.XYmax);
                else
                    fprintf('\nWarning: psi has been initialized.\n\n');
                end
                
                %initialize record psi time two_component
                if psiobj.record_psi~=0
                    if  isempty(psiobj.TotalTime)
                        psiobj.re_psi_end=psiobj.TotalTime;
                        fprintf('\nWarning: TotalTime is not defined.\n\n');
                    else
                        psiobj.record_psi_t=psiobj.re_psi_begin:psiobj.re_psi_dt:psiobj.re_psi_end;
                        %                         if psiobj.component_number==1
                        %                             psiobj.psi_trim=zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t));
                        %                         elseif psiobj.component_number==2;
                        psiobj.psiL_trim=complex(zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t)),0);
                        psiobj.psiR_trim=complex(zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t),0));
                        psiobj.nL_trim=zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t));
                        psiobj.nR_trim=zeros(psiobj.N,psiobj.N,length(psiobj.record_psi_t));
                        %                         end
                        psiobj.TotalStep=floor(psiobj.TotalTime/psiobj.dt);
                    end
                end
                
            end
            
            %initialize Time and Energy (two and single component)
            if isempty(psiobj.Time) && isempty(psiobj.Energy)
                psiobj.E_step_interval=floor(psiobj.E_t_interval/psiobj.dt);
                psiobj.Time=zeros(1,floor(psiobj.TotalStep/psiobj.E_step_interval));
                psiobj.Energy=zeros(1,floor(psiobj.TotalStep/psiobj.E_step_interval));
                psiobj.Energy_unmal=psiobj.Energy;
                psiobj.Psi_int=psiobj.Energy;
                psiobj.Lz=psiobj.Energy;
                psiobj.mod_psiLR_square=zeros(1,floor(psiobj.TotalStep/psiobj.E_step_interval));
                psiobj.sx_int=psiobj.Energy;
                psiobj.sy_int=psiobj.Energy;
                psiobj.sz_int=psiobj.Energy;
                psiobj.step_E=1;
            else
                fprintf('\nError: Time or Energy has been initialized.\n\n');
            end
            
            %initialize integration area
            psiobj.int_area=heaviside1(psiobj.int_r-sqrt(psiobj.x.^2+psiobj.y.^2));
            
        end %initialization function end
        
        function record_Psi(psiobj)
            if psiobj.component_number==1
               if psiobj.record_psi==0
                    fprintf('\nWarning: will not be recorded.\n\n');
                end
                
                psiobj.psi_trim(:,:,psiobj.step_record_psi)=psiobj.psi;
                psiobj.n_trim(:,:,psiobj.step_record_psi)=psiobj.n;
                fprintf('               psi recorded\n');
                psiobj.step_record_psi=psiobj.step_record_psi+1;
                
            elseif psiobj.component_number==2
                if psiobj.record_psi==0
                    fprintf('\nWarning: will not be recorded.\n\n');
                end
                
                psiobj.psiL_trim(:,:,psiobj.step_record_psi)=psiobj.psiL;
                psiobj.psiR_trim(:,:,psiobj.step_record_psi)=psiobj.psiR;
                psiobj.nL_trim(:,:,psiobj.step_record_psi)=psiobj.nL;
                psiobj.nR_trim(:,:,psiobj.step_record_psi)=psiobj.nR;
                
                fprintf('               psi recorded\n');
                psiobj.step_record_psi=psiobj.step_record_psi+1;
            end
        end %record psi
            
        function record_Time_Energy(psiobj)
            %record Time
            psiobj.Time(psiobj.step_E)=double(psiobj.dt*psiobj.step_i);
            %calculate the energy term
            
            if psiobj.component_number==1 %Single component
                [psix, psiy]=gradient(psiobj.psi,psiobj.hspace);
                mod_psi_square=abs(psiobj.psi).^2;
                
                %store psi_int
                psiobj.Psi_int(psiobj.step_E)=(psiobj.hspace^2*sum(sum(mod_psi_square.*psiobj.int_area)));
                
                %store Lz (normalized)
                psiobj.Lz(psiobj.step_E)=psiobj.hspace^2*sum(sum((-1i)*conj(psiobj.psi).*(psiobj.x.*psiy-psiobj.y.*psix).*psiobj.int_area))/psiobj.Psi_int(psiobj.step_E);
                
                %calculate Energy_xy; real part of Eq.(2.10) in Bao's paper
                Energy_xy=0.5*(abs(psix).^2+abs(psiy).^2)+psiobj.gR*psiobj.n.*mod_psi_square+(0.5)*psiobj.ua*mod_psi_square.^2;
                
                %store Energy term (un-normalized)
                psiobj.Energy_unmal(psiobj.step_E)=(psiobj.hspace^2*sum(sum(Energy_xy.*psiobj.int_area)));

                %store Energy term (normalized)
                psiobj.Energy(psiobj.step_E)=psiobj.Energy_unmal(psiobj.step_E)/psiobj.Psi_int(psiobj.step_E);
                
            elseif psiobj.component_number==2 %two component
                hspace=psiobj.hspace;
                [psiL_x, psiL_y]=gradient(psiobj.psiL,psiobj.hspace);
                [psiR_x, psiR_y]=gradient(psiobj.psiR,psiobj.hspace);
                mod_psiL=abs(psiobj.psiL);
                mod_psiR=abs(psiobj.psiR);
                mod_psiL_square=mod_psiL.^2;
                mod_psiR_square=mod_psiR.^2;
                mod_psiLR_square_temp=mod_psiL+mod_psiR_square;
                
                %update integration area
                if psiobj.fixed_int_r==0 %only for round Gaussin pump and pulse
                    for it=psiobj.N/2:psiobj.N
                        if mod_psiLR_square_temp(psiobj.N/2+1,it)<=psiobj.int_tolerant
                           psiobj.int_r=(it-psiobj.N/2)*psiobj.hspace;
                           break;
                        end
                        if it==psiobj.N
                           psiobj.int_r=psiobj.hspace;
                        end
                    end        
                    psiobj.int_area=heaviside1(psiobj.int_r-sqrt(psiobj.x.^2+psiobj.y.^2));
                end
                
                if psiobj.Dk==0 %no TE-TM
                    Energy_xy=0.5*(abs(psiL_x).^2+abs(psiL_y).^2+abs(psiR_x).^2+abs(psiR_y).^2)+psiobj.gR*psiobj.nL.*mod_psiL_square+psiobj.gR*psiobj.nR.*mod_psiR_square +0.5*psiobj.ua*(mod_psiL.^4+mod_psiR.^4)+psiobj.ub*mod_psiL_square.*mod_psiR_square+2*psiobj.J*real(psiobj.psiL.*conj(psiobj.psiR));
                else %TE-TM
                    [psiL_xx,psiL_xy]=gradient(psiL_x,hspace);
                    [psiR_xx,psiR_xy]=gradient(psiR_x,hspace);
                    [~,psiL_yy]=gradient(psiL_y,hspace);
                    [~,psiR_yy]=gradient(psiR_y,hspace);
                    Energy_xy=0.5*(abs(psiL_x).^2+abs(psiL_y).^2+abs(psiR_x).^2+abs(psiR_y).^2)+psiobj.gR*psiobj.nL.*mod_psiL_square+psiobj.gR*psiobj.nR.*mod_psiR_square +psiobj.ua*(mod_psiL.^4+mod_psiR.^4)+2*psiobj.ub*mod_psiL_square.*mod_psiR_square+2*psiobj.J*real(psiobj.psiL.*conj(psiobj.psiR))...
                              +psiobj.Dk*(conj(psiobj.psiR).*(-psiL_xx+2i*psiL_xy+psiL_yy)+conj(psiobj.psiL).*(-psiR_xx-2i*psiR_xy+psiR_yy));
                    if mod(psiobj.step_E,20)==0
                       fprintf('TE-TM energy recorded.\n')
                    end
                end
                
                %record Energy (integrated)
                psiobj.Energy(psiobj.step_E)=psiobj.hspace^2*sum(sum(Energy_xy.*psiobj.int_area));
                %record psiLR square (integrated)
                psiobj.mod_psiLR_square(psiobj.step_E)=psiobj.hspace^2*sum(sum(mod_psiLR_square_temp.*psiobj.int_area));
                
            end
            
            %record integrated Strokes parameters
            if psiobj.component_number==2 && psiobj.record_sxyz_int~=0
                theta_temp=angle(psiobj.psiL)-angle(psiobj.psiR);
                
                Sx=mod_psiR.*mod_psiL.*cos(theta_temp);
                Sy=mod_psiR.*mod_psiL.*sin(theta_temp);
                Sz=0.5*(mod_psiR.^2-mod_psiL.^2);
                
                Sx_int=psiobj.hspace^2*sum(sum(Sx.*psiobj.int_area));
                Sy_int=psiobj.hspace^2*sum(sum(Sy.*psiobj.int_area));
                Sz_int=psiobj.hspace^2*sum(sum(Sz.*psiobj.int_area));
                
                S0=sqrt(Sx_int^2+Sy_int^2+Sz_int^2);
                
                psiobj.sx_int(psiobj.step_E)=Sx_int/S0;
                psiobj.sy_int(psiobj.step_E)=Sy_int/S0;
                psiobj.sz_int(psiobj.step_E)=Sz_int/S0;

            end
            
            psiobj.step_E=psiobj.step_E+1;
            
        end %record Time_Energy function end
        
        function clear_all(psiobj)
            if psiobj.component_number==1
                psiobj.psi=[];psiobj.n=[];
            elseif psiobj.component_number==2
                psiobj.nL=[]; psiobj.nR=[];
                psiobj.psiL=[];psiobj.psiR=[];
            end
            psiobj.x=[]; psiobj.y=[];
            psiobj.Time=[];psiobj.Energy=[];
            psiobj.mod_psiLR_square=[];
            psiobj.sx_int=[];
            psiobj.sy_int=[];
            psiobj.sz_int=[];
            psiobj.step_E=1; psiobj.step_i=[];
        end
        
    end
    
end









