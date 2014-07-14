function     Gaussian_O_FollowUp_Co_HalfVortex_Function(Pmp,Pse,Psi,M)

%Full simulation program with Coherent pulse (Object oriented)
%save Time Energy
%record mov
%record trimmed psi
fprintf('\nGaussian FollowUp program.\n')
tic;

%set dataFolder name
currentFolder = pwd;
folder_name_parameter=sprintf('IR=%.1f_J=%.1f',Psi.IR,Psi.J);
dataFolder=strcat('GauPsi_',M.Simulation_version,'_CoPulse_',folder_name_parameter);

load(strcat(currentFolder,'\',dataFolder,'\','Gau_Standard_Steady_State_',M.Simulation_version,'_m0_',num2str(Psi.Ini_m),'.mat'));
%Check parameters insistent
if ~Pmp.Check_Parameter(Ori_Pump_Parameter);
   fprintf('\nError: Parameters do not fit.\n\n');
   return;
end

N=Pmp.N;
XYmax=Pmp.XYmax;
r_damping=XYmax-14; %radius for top-hat damping
a_damping=0.3; %slope for the damping mesa

TotalTime=Pmp.TotalTime;
dt=Pmp.dt;
if TotalTime<=max(Time);
   fprintf('\nError: TotalTime is smaller then steady state time\n\n');
   return;
end

%Coefficients in the ODGPE: (parameters from the PRA)
ua=Pmp.ua; %same spin polariton-polariton scattering
ub=Pmp.ub; %cross spin scattering
gR=Pmp.gR; %LP interaction with incoherent reservoir
GammaC=Pmp.GammaC;
GammaR=Pmp.GammaR; % GammaR_ratio=GammaR/GammaC
R=Pmp.R; %stimulated scattering rate
J=Pmp.J; %general coefficient
IR=Psi.IR;

%initialize psi
Psi.initialize();
psiL=Finial_State(:,:,1);
psiR=Finial_State(:,:,2);
nL=Finial_State(:,:,3);
nR=Finial_State(:,:,4);
Psi.step_E=length(Time);
Psi.Time(1,1:Psi.step_E)=Time;
Psi.Energy(1,1:Psi.step_E)=Energy;
Psi.mod_psiLR_square(1,1:Psi.step_E)=mod_psiLR_square;
if Psi.record_sxyz_int~=0
    if ~isempty(sx_int)
        Psi.sx_int(1,1:Psi.step_E)=sx_int;
        Psi.sy_int(1,1:Psi.step_E)=sy_int;
        Psi.sz_int(1,1:Psi.step_E)=sz_int;
    end
end
step_0=floor(max(Time)/dt);

%construct the loaded Pump
Pmp.Build_u();
PumpL=Pmp.PL*abs(Pmp.u).^2;
PumpR=(1-Pmp.PL)*abs(Pmp.u).^2;

%construct Pulse
Pse.Build_u();
Pulse_add_flag=0; %0:pump only, 1:pump+pulse, 2:pulse ended

%construct lattice size and k vector space
x=Pmp.x;
y=Pmp.y;
kvector=fftshift(-N/2:N/2-1)*2*pi/(2*XYmax);
[kx, ky]=meshgrid(kvector, kvector);
kinetic_factor_quarter=exp((-1i*(kx.^2+ky.^2)*dt)/4);

%Absorbing boundaries (damping)
damping=0.25*(1+tanh(a_damping*((x.^2+y.^2).^0.5+r_damping))).*(1+tanh(a_damping*(-(x.^2+y.^2).^0.5+r_damping)));

%Saved data file name
file_name_parameter=sprintf('Pmp0=%.2f_PmpPL=%.1f_PmpSmax=%d_PmpSmay=%d_N=%d_T=%d',Pmp.Pbar0,Pmp.PL,Pmp.Sigma_x,Pmp.Sigma_y,N,floor(TotalTime));
Mov_file_name_psi=strcat('GauPsi_',M.Simulation_version,'_CoPulse_',file_name_parameter);
Mov_file_name_sxyz=strcat('Gausxyz_',M.Simulation_version,'_CoPulse_',file_name_parameter);
file_name=strcat('GauPsi_',M.Simulation_version,'_Co_',file_name_parameter);

if Pse.force_pulse_time<=TotalTime
    file_name_parameter=sprintf('Pmp0=%.1f_PmpPL=%.1f_PmpSx=%d_PmpSy=%d_Pse0=%.2f_PseP=%d_PsePL=%.2f_PseSx=%d_PseSy=%d_Px0=%d_Py0=%d_N=%d_T=%d',...
                                Pmp.Pbar0,Pmp.PL,Pmp.Sigma_x,Pmp.Sigma_y,Pse.Pbar0,round(pi/Pse.C_phase),Pse.PL,Pse.Sigma_x,Pse.Sigma_y,Pse.x0,Pse.y0,N,floor(TotalTime));
    Mov_file_name_psi=strcat('GauPsi_',M.Simulation_version,'_CoPulse_',file_name_parameter);
    Mov_file_name_sxyz=strcat('Gausxyz_',M.Simulation_version,'_CoPulse_',file_name_parameter);
    file_name=strcat('GauPsi_',M.Simulation_version,'_Co_',file_name_parameter);
end

%GPU calculation initializing
if M.GPU_calculation~=0
    GPU_Device = gpuDevice;
    gpuDevice(1);
    kinetic_factor_quarter_gpu=gpuArray(kinetic_factor_quarter);
    damping_gpu=gpuArray(damping);
    if M.GPU_calculation==2 %full GPU calculation
        psiL_gpu=gpuArray(psiL);
        psiR_gpu=gpuArray(psiR);
        nL_gpu=gpuArray(nL);
        nR_gpu=gpuArray(nR);
        PumpL_gpu=gpuArray(PumpL);
        PumpR_gpu=gpuArray(PumpR);
    end
    fprintf('\nGPU calculation enabled\n');
end

fftw('planner','exhaustive');

fprintf('\nCalculation begins\n\n');

for step_i=step_0:floor(TotalTime/dt)
   %Core calculation----------------------------------
    if M.GPU_calculation==2 %Full GPU calculation
        %reservoirs update
        nL_gpu=nL_gpu.*exp(-(GammaR+R*abs(psiL_gpu).^2)*dt)+PumpL_gpu*dt;
        nR_gpu=nR_gpu.*exp(-(GammaR+R*abs(psiR_gpu).^2)*dt)+PumpR_gpu*dt;
        
        psiL1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiL_gpu)));
        psiR1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiR_gpu)));
        
        if Pulse_add_flag~=1
            psiL_temp_gpu=psiL1_gpu;
            psiL1_gpu=psiL1_gpu.*exp((ua*abs(psiL1_gpu).^2+ub*abs(psiR1_gpu).^2+gR*(nL_gpu)+0.5i*(R*nL_gpu-GammaC))*(-1i)*dt)-1i*J*psiR1_gpu*dt;
            psiR1_gpu=psiR1_gpu.*exp((ua*abs(psiR1_gpu).^2+ub*abs(psiL_temp_gpu).^2+gR*(nR_gpu)+0.5i*(R*nR_gpu-GammaC))*(-1i)*dt)-1i*J*psiL_temp_gpu*dt;
        else
            psiL_temp_gpu=psiL1_gpu; %pulse added
            psiL1_gpu=psiL1_gpu.*exp((ua*abs(psiL1_gpu).^2+ub*abs(psiR1_gpu).^2+gR*(nL_gpu)+0.5i*(R*nL_gpu-GammaC))*(-1i)*dt)-1i*J*psiR1_gpu*dt-1i*gpuArray(Pse.u)*dt;
            psiR1_gpu=psiR1_gpu.*exp((ua*abs(psiR1_gpu).^2+ub*abs(psiL_temp_gpu).^2+gR*(nR_gpu)+0.5i*(R*nR_gpu-GammaC))*(-1i)*dt)-1i*J*psiL_temp_gpu*dt-1i*gpuArray(Pse.u)*dt;
        end
        
        psiL_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiL1_gpu))).*damping_gpu;
        psiR_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiR1_gpu))).*damping_gpu;
        
    else %CPU or Half GPU
        %reservoirs update
        nL=nL.*exp(-(GammaR+R*abs(psiL).^2)*dt)+PumpL*dt;
        nR=nR.*exp(-(GammaR+R*abs(psiR).^2)*dt)+PumpR*dt;
        
        %advance spacially the first time
        if M.GPU_calculation==1 %GPU FFT
            psiL1=gather(ifft2(kinetic_factor_quarter_gpu.*fft2(gpuArray(psiL))));
            psiR1=gather(ifft2(kinetic_factor_quarter_gpu.*fft2(gpuArray(psiR))));
        else
            psiL1=ifft2((kinetic_factor_quarter.*fft2(psiL)));
            psiR1=ifft2((kinetic_factor_quarter.*fft2(psiR)));
        end
        
        %update psiL1,R1 by the pseudo-potential
        if Pulse_add_flag~=1
            psiL_temp=psiL1;
            psiL1=psiL1.*exp((ua*abs(psiL1).^2+ub*abs(psiR1).^2+gR*(nL)+0.5i*(R*nL-GammaC))*(-1i)*dt)-1i*J*psiR1*dt;
            psiR1=psiR1.*exp((ua*abs(psiR1).^2+ub*abs(psiL_temp).^2+gR*(nR)+0.5i*(R*nR-GammaC))*(-1i)*dt)-1i*J*psiL_temp*dt;
        else
            psiL_temp=psiL1;
            psiL1=psiL1.*exp((ua*abs(psiL1).^2+ub*abs(psiR1).^2+gR*(nL)+0.5i*(R*nL-GammaC))*(-1i)*dt)-1i*J*psiR1*dt-1i*Pse.u*dt;
            psiR1=psiR1.*exp((ua*abs(psiR1).^2+ub*abs(psiL_temp).^2+gR*(nR)+0.5i*(R*nR-GammaC))*(-1i)*dt)-1i*J*psiL_temp*dt-1i*Pse.u*dt;
        end
        
        %advance the second time with boundary condition (damping)
        if M.GPU_calculation==1 %GPU FFT
            psiL=gather(ifft2((kinetic_factor_quarter_gpu.*fft2(gpuArray(psiL1)))).*damping_gpu);
            psiR=gather(ifft2((kinetic_factor_quarter_gpu.*fft2(gpuArray(psiR1)))).*damping_gpu);
        else
            psiL=ifft2((kinetic_factor_quarter.*fft2(psiL1))).*damping;
            psiR=ifft2((kinetic_factor_quarter.*fft2(psiR1))).*damping;
        end
    end
    %----------------------------------

    time=double(step_i*dt);
    
    if Pulse_add_flag==1
        fprintf('t = %.3f      Pulse added\n',time);
    else   fprintf('t = %.3f\n',time);
    end
    
    if M.GPU_calculation==2
        if (mod(step_i,Psi.E_step_interval)==0)...
                || (Psi.record_psi~=0 && time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0))...
                || (M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0))
            psiL=gather(psiL_gpu);
            psiR=gather(psiR_gpu);
        end
    end
    
    %record Time Energy, mod_psiLR, and sxyz_int
    if mod(step_i,Psi.E_step_interval)==0
        Psi.step_i=step_i;
        Psi.psiL=psiL;
        Psi.psiR=psiR;
        Psi.record_Time_Energy();
    end
    
    %store trimmed psi
    if Psi.record_psi~=0
        if time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0)
            Psi.psiL=psiL;
            Psi.psiR=psiR;
            Psi.record_Psi();
        end
    end
    
    %record Psi Mov
    if M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0)
        
        psi_h=(psiR+psiL)/sqrt(2);
        psi_v=1i*(psiR-psiL)/sqrt(2);
        
        mod_psi_h_square=abs(psi_h).^2;
        mod_psi_v_square=abs(psi_v).^2;
        phi_h=wrapToPi(angle(psi_h));
        phi_v=wrapToPi(angle(psi_v));
        
        mod_psiL_square=abs(psiL).^2;
        mod_psiR_square=abs(psiR.^2);
        phiL=wrapToPi(angle(psiL));
        phiR=wrapToPi(angle(psiR));
        
        if M.record_mov_psi~=0
            
            if M.re_psi_initialized==0 %psi windows and mesh initialization
                
                psi_fig=figure('position',[50 100 1860 700],'renderer','zbuffer');
                M.re_psi_initialized=1;
                M.mov_file_name_psi=Mov_file_name_psi;
                M.currentFolder=currentFolder;
                M.dataFolder=dataFolder;
                
                psi_h_plot=subplot(2,4,1);
                mesh_psi_h=mesh(x,y,mod_psi_h_square);
                axis square;
                set(mesh_psi_h,'ZDataSource','mod_psi_h_square');
                set(psi_h_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                if Pulse_add_flag==1%Pse.force_pulse_time>=time && (Pse.force_pulse_time+Pse.pulse_totaltime)<=time
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_h_plot,{st;'|\psi_h|^2'});
                
                psi_v_plot=subplot(2,4,2);
                mesh_psi_v=mesh(x,y,mod_psi_v_square);
                axis square;
                set(mesh_psi_v,'ZDataSource','mod_psi_v_square');
                set(psi_v_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('|\psi_v|^2');
                
                               
                %plot phase
                phi_h_plot=subplot(2,4,5);
                mesh_phi_h=mesh(x,y,phi_h);
                axis square;
                set(mesh_phi_h,'ZDataSource','phi_h');
                set(phi_h_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_h');
                
                phi_v_plot=subplot(2,4,6);
                mesh_phi_v=mesh(x,y,phi_v);
                axis square;
                set(mesh_phi_v,'ZDataSource','phi_v');
                set(phi_v_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_v');
                
                subplot(2,4,3);
                mesh_psiL=mesh(x,y,mod_psiL_square);
                axis square;
                set(mesh_psiL,'ZDataSource','mod_psiL_square');
                set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('|\psi_L|^2');
                
                subplot(2,4,4);
                mesh_psiR=mesh(x,y,mod_psiR_square);
                axis square;
                set(mesh_psiR,'ZDataSource','mod_psiR_square');
                set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('|\psi_R|^2');
                
                subplot(2,4,7);
                mesh_phiL=mesh(x,y,phiL);
                axis square;
                set(mesh_phiL,'ZDataSource','phiL');
                set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_L');
                
                subplot(2,4,8);
                mesh_phiR=mesh(x,y,phiR);
                axis square;
                set(mesh_phiR,'ZDataSource','phiR');
                set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_R');
                
                
            else
                refreshdata(mesh_psi_h,'caller');
                refreshdata(mesh_psi_v,'caller');
                refreshdata(mesh_psiL,'caller');
                refreshdata(mesh_psiR,'caller');
                refreshdata(mesh_phi_h,'caller');
                refreshdata(mesh_phi_v,'caller');
                refreshdata(mesh_phiL,'caller');
                refreshdata(mesh_phiR,'caller');
                
                if Pulse_add_flag==1
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_h_plot,{st;'|\psi_h|^2'});
            end
            M.Mov_psi(psi_fig);
        end % psi movie end
        
        if M.record_mov_sxyz~=0 %record sxyz movie
            
            %calculate Strokes parameters
            if isempty(mod_psiL)
                mod_psiL=abs(psiL);
                mod_psiR=abs(psiR);
                phiL=wrapToPi(angle(psiL));
                phiR=wrapToPi(angle(psiR));
                theta=phiL-phiR;
                mod_psiLR_square=mod_psiL.^2+mod_psiR.^2;
            end
            
            nc=mod_psiLR_square;
            sx=mod_psiR.*mod_psiL.*cos(theta)./nc;
            sy=mod_psiR.*mod_psiL.*sin(theta)./nc;
            sz=0.5*(mod_psiR.^2-mod_psiL.^2)./nc;
            
            for i_n=1:N %non-number error check
                for j_n=1:N
                    if isnan(sx(i_n,j_n))
                        sx(i_n,j_n)=0;
                    end
                    if isnan(sy(i_n,j_n))
                        sy(i_n,j_n)=0;
                    end
                    if isnan(sz(i_n,j_n))
                        sz(i_n,j_n)=0;
                    end
                end
            end
            
            if M.re_sxyz_initialized==0
                sxyz_fig=figure('position',[50 100 1100 900],'renderer','zbuffer');
                M.mov_file_name_sxyz=Mov_file_name_sxyz;
                M.currentFolder=currentFolder;
                M.dataFolder=dataFolder;
                M.re_sxyz_initialized=1;
                
                axis_sx=subplot(2,2,1);
                mesh_sx=mesh(x,y,sx);
                set(mesh_sx,'ZDataSource','sx');
                axis square;
                set(axis_sx,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-0.5 0.5],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                caxis([-0.5 0.5]);
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar
                if Pulse_add_flag==1
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(axis_sx,{st;'s_x'});
                
                axis_sy=subplot(2,2,2);
                mesh_sy=mesh(x,y,sy);
                set(mesh_sy,'ZDataSource','sy');
                axis square;
                set(axis_sy,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-0.5 0.5],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                caxis([-0.5 0.5]);
                title('s_y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar
                
                axis_sz=subplot(2,2,3);
                mesh_sz=mesh(x,y,sz);
                set(mesh_sz,'ZDataSource','sz');
                axis square;
                set(axis_sz,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-0.5 0.5],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                title('s_z');
                caxis([-0.5 0.5]);
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar
                
                % I_J mesh
                if J~=0
                    I_J_plot=subplot(2,2,4);
                    mesh_I_J=mesh(x,y,I_J);
                    axis square;
                    set(mesh_I_J,'ZDataSource','I_J');
                    set(I_J_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('I_J');
                end
                
            else
                refreshdata(mesh_sx,'caller');
                refreshdata(mesh_sy,'caller');
                refreshdata(mesh_sz,'caller');
                
                if J~=0
                    refreshdata(mesh_I_J,'caller');
                end
                
                if Pulse_add_flag==1
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(axis_sx,{st;'s_x'});
            end
            
            M.Mov_sxyz(sxyz_fig);
            
        end %record sxyz mov end
    end %record two movies end
    
    %Coherent pulse
    if Pulse_add_flag==0 %stage 1
        if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
            fprintf('               pulse begins\n');
            Pulse_add_flag=1;%enter state 2
        end
    elseif Pulse_add_flag==1 %stage 2
            
            %Pse.time=double(time-Pse.force_pulse_time);
            %adiabatic_term=0.5*(1+tanh(Pse.ada_coe*(Pse.time-Pse.ada_shift)));
            
        if (Pse.force_pulse_time+Pse.pulse_totaltime)<=time %stage 3
            fprintf('               pulse ends\n')
            Pse.clear_all();
            Pulse_add_flag=2;%enter state 3
        end
    end %pulse update end
    
end %end of main loop

%clean up movie writer
if M.record_movie~=0
    if ~isempty(M.vidObjpsi)
        close(M.vidObjpsi);
        close(psi_fig);
    end
    if ~isempty(M.vidObjsxyz)
        close(M.vidObjsxyz);
        close(sxyz_fig);
    end
end
%clean up GPU calculation parameters
if M.GPU_calculation~=0
    if M.GPU_calculation==2
       psiL=gather(psiL_gpu);
       psiR=gather(psiR_gpu);
       nL=gather(nL_gpu);
       nR=gather(nR_gpu);
    end
   delete(GPU_Device);    
end

%draw final state and Energy evolution
fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','zbuffer','Visible','off');

psi_h=(psiR+psiL)/sqrt(2);
psi_v=1i*(psiR-psiL)/sqrt(2);

subplot(2,4,1);
mesh(x,y,abs(psi_h).^2);
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.6427;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title({sprintf('               t=%.1f',step_i*dt);'|\psi_h|^2'});

subplot(2,4,2);
mesh(x,y,abs(psi_v).^2);
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.6427;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi_v|^2');

subplot(2,4,3);
mesh(x,y,abs(psiL).^2);
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.6427;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi_L|^2');

subplot(2,4,4);
mesh(x,y,abs(psiR).^2);
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.6427;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi_R|^2');

%phase
subplot(2,4,5);
mesh(x,y,wrapToPi(angle(psi_h)));
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.32;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
caxis([-pi pi]);
colorbar;
title('\phi_h');

a_phi_v=subplot(2,4,6);
mesh(x,y,wrapToPi(angle(psi_v)));
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.32;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
caxis([-pi pi]);
colorbar;
title('\phi_v');

a_phiL=subplot(2,4,7);
mesh(x,y,wrapToPi(angle(psiL)));
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.32;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
caxis([-pi pi]);
colorbar;
title('\phi_L');

subplot(2,4,8);
mesh(x,y,wrapToPi(angle(psiR)));
axis square;
set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
po_temp=get(gca,'position');
po_temp(2)=0.32;
set(gca,'position',po_temp);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
caxis([-pi pi]);
colorbar;
title('\phi_R');

a_E=axes('position',[0.0961    0.0805    0.1800 0.23]);
plot(Psi.Time,Psi.Energy);
axis square;
xlabel('t');
ylabel('Energy');
set(gca,'fontsize',Psi.font_size);

saveas(fig_E,strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'),'png');
%print(fig_E,'-zbuffer','-dpng','-r250',strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'));
close(fig_E);

fprintf('Storing Data\n ');

Time=Psi.Time;
Energy=Psi.Energy;
record_psi_t=Psi.record_psi_t;
psiL_trim=Psi.psiL_trim;
psiR_trim=Psi.psiR_trim;
Pmp.clear_all();
Pse.clear_all();
Psi.clear_all();
M.clear_all();
Ori_Pump_Parameter=Pmp;
save(strcat(currentFolder,'\',dataFolder,'\',file_name,'.mat'),'Ori_Pump_Parameter','record_psi_t','psiL_trim','psiR_trim','Time','Energy','-mat','-v7.3');

fprintf('\nSimulation function finished.\n\n')
toc;

end



