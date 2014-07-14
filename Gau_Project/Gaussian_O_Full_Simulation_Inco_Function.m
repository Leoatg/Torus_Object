function     Gaussian_O_Full_Simulation_Inco_Function(Pmp,Pse,Psi,M)

%Full simulation program with Incoherent pulse (Object oriented)
%save Time Energy
%record mov
%record trimmed psi
fprintf('\nFull simulation program.\n')
tic;
N=Pmp.N;
XYmax=Pmp.XYmax;

TotalTime=Pmp.TotalTime;
dt=Pmp.dt;

%Coefficients in the ODGPE: (parameters from the PRA)
ua=Pmp.ua; %same spin polariton-polariton scattering
ub=Pmp.ub; %cross spin scattering
gR=Pmp.gR; %LP interaction with incoherent reservoir
GammaC=Pmp.GammaC;
GammaR=Pmp.GammaR; % GammaR_ratio=GammaR/GammaC
R=Pmp.R; %stimulated scattering rate
J=Pmp.J; %general coefficient
IR=Psi.IR;

%construct Pump
Pmp.Build_u();
PumpL=Pmp.PL*abs(Pmp.u).^2;
PumpR=(1-Pmp.PL)*abs(Pmp.u).^2;

%construct Pulse
Pse.Build_u();
Pulse_add_flag=0; %0:pump only, 1:pump+pulse, 2:pulse ended

%initialize psi
Psi.initialize();
psiL=Psi.psiL;
psiR=Psi.psiR;
nL=Psi.nL;
nR=Psi.nR;

%construct lattice size and k vector space
x=Pmp.x;
y=Pmp.y;
kvector=fftshift(-N/2:N/2-1)*2*pi/(2*XYmax);
[kx, ky]=meshgrid(kvector, kvector);
kinetic_factor_quarter=exp((-1i*(kx.^2+ky.^2)*dt)/4);

%Absorbing boundaries (damping)
r_damping=XYmax-Psi.damping_shift;%XYmax-14; %radius for top-hat damping
a_damping=Psi.a_damping; %slope for the damping mesa
damping=0.25*(1+tanh(a_damping*((x.^2+y.^2).^0.5+r_damping))).*(1+tanh(a_damping*(-(x.^2+y.^2).^0.5+r_damping)));
x_axis=-XYmax:Pmp.hspace:XYmax;

%Saved data file name
file_name_parameter=sprintf('Pmp0=%.2f_PmpPL=%.3f_PmpSmax=%d_PmpSmay=%d_XYmax=%d_N=%d_T=%d',Pmp.Pbar0,Pmp.PL,Pmp.Sigma_x,Pmp.Sigma_y,Pmp.XYmax,N,floor(TotalTime));
Mov_file_name_psi=strcat('GauPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
Mov_file_name_sxyz=strcat('Gausxyz_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
file_name=strcat('GauPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);

if Pse.force_pulse_time<=TotalTime
    file_name_parameter=sprintf('Pmp0=%.2f_PmpPL=%.2f_PmpSmax=%d_PmpSmay=%d_Pse0=%.2f_PsePL=%.2f_PseSmax=%d_PseSmay=%d_N=%d_T=%d',...
                                Pmp.Pbar0,Pmp.PL,Pmp.Sigma_x,Pmp.Sigma_y,Pse.Pbar0,Pse.PL,Pse.Sigma_x,Pse.Sigma_y,N,floor(TotalTime));
    Mov_file_name_psi=strcat('GauPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
    Mov_file_name_sxyz=strcat('Gausxyz_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
    file_name=strcat('GauPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
end

currentFolder = pwd;
folder_name_parameter=sprintf('IR=%.1f_J=%.1f',IR,J);
dataFolder=strcat('GauPsi_',M.Simulation_version,'_IncoPulse_',folder_name_parameter);

if exist(dataFolder,'dir')==0
    mkdir(dataFolder);
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

for step_i=1:floor(TotalTime/dt)
    %Core calculation----------------------------------
    if M.GPU_calculation==2 %Full GPU calculation
        %reservoirs update
        nL_gpu=nL_gpu.*exp(-(GammaR+R*abs(psiL_gpu).^2)*dt)+PumpL_gpu*dt;
        nR_gpu=nR_gpu.*exp(-(GammaR+R*abs(psiR_gpu).^2)*dt)+PumpR_gpu*dt;
        
        psiL1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiL_gpu)));
        psiR1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psiR_gpu)));
        
        psiL_temp_gpu=psiL1_gpu;
        psiL1_gpu=psiL1_gpu.*exp((ua*abs(psiL1_gpu).^2+ub*abs(psiR1_gpu).^2+gR*(nL_gpu)+0.5i*(R*nL_gpu-GammaC))*(-1i)*dt)-1i*J*psiR1_gpu*dt;
        psiR1_gpu=psiR1_gpu.*exp((ua*abs(psiR1_gpu).^2+ub*abs(psiL_temp_gpu).^2+gR*(nR_gpu)+0.5i*(R*nR_gpu-GammaC))*(-1i)*dt)-1i*J*psiL_temp_gpu*dt;
        
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
        psiL_temp=psiL1;
        psiL1=psiL1.*exp((ua*abs(psiL1).^2+ub*abs(psiR1).^2+gR*(nL)+0.5i*(R*nL-GammaC))*(-1i)*dt)-1i*J*psiR1*dt;
        psiR1=psiR1.*exp((ua*abs(psiR1).^2+ub*abs(psiL_temp).^2+gR*(nR)+0.5i*(R*nR-GammaC))*(-1i)*dt)-1i*J*psiL_temp*dt;
        
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
        fprintf('t = %.5f      Pulse added\n',time);
    else   fprintf('t = %.5f\n',time);
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
        if M.GPU_calculation==2
            Psi.nL=gather(nL_gpu);
            Psi.nR=gather(nR_gpu);
        else 
            Psi.nL=nL;
            Psi.nR=nR;
        end
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
        
        mod_psiL=abs(psiL);
        mod_psiR=abs(psiR);
        mod_psiL_square=mod_psiL.^2;
        mod_psiR_square=mod_psiR.^2;
        mod_psiLR_square=mod_psiL_square+mod_psiR_square;
        phiL=wrapToPi(angle(psiL));
        phiR=wrapToPi(angle(psiR));
        theta=phiL-phiR;
        I_J=2*J*mod_psiL.*mod_psiR.*sin(theta);
        
        if M.record_mov_psi~=0
            
            if M.re_psi_initialized==0 %psi windows and mesh initialization
                
                psi_fig=figure('position',[50 80 1800 900],'renderer','zbuffer');
                M.re_psi_initialized=1;
                M.mov_file_name_psi=Mov_file_name_psi;
                M.currentFolder=currentFolder;
                M.dataFolder=dataFolder;
                
                psi_L_plot=subplot(2,3,1);
                mesh_psiL=mesh(x,y,mod_psiL_square);
                axis square;
                set(mesh_psiL,'ZDataSource','mod_psiL_square');
                set(psi_L_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                if Pulse_add_flag==1%Pse.force_pulse_time>=time && (Pse.force_pulse_time+Pse.pulse_totaltime)<=time
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_L_plot,{st;'|\psi_L|^2'});
                
                psi_R_plot=subplot(2,3,2);
                mesh_psiR=mesh(x,y,mod_psiR_square);
                axis square;
                set(mesh_psiR,'ZDataSource','mod_psiR_square');
                set(psi_R_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('|\psi_R|^2');
                
                psi_LR_plot=subplot(2,3,3);
                mesh_psi_LR=mesh(x,y,mod_psiLR_square);
                axis square;
                set(mesh_psi_LR,'ZDataSource','mod_psiLR_square');
                set(psi_LR_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('|\psi_L|^2 + |\psi_R|^2');
                
                %plot phase
                phi_L_plot=subplot(2,3,4);
                mesh_phi_L=mesh(x,y,phiL);
                axis square;
                set(mesh_phi_L,'ZDataSource','phiL');
                set(phi_L_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_L');
                
                phi_R_plot=subplot(2,3,5);
                mesh_phi_R=mesh(x,y,phiR);
                axis square;
                set(mesh_phi_R,'ZDataSource','phiR');
                set(phi_R_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                caxis([-pi pi]);
                colorbar;
                title('\phi_R');
                
%                 axis_pump=subplot(2,3,6);
%                 mesh_pumpL=mesh(axis_pump,x,y,PumpL);
%                 set(mesh_pumpL,'ZDataSource','PumpL');
%                 axis square;
%                 set(axis_pump,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
%                 view(0,90);
%                 zoom(Psi.zoom_factor);
%                 colorbar
%                 xlabel('x');
%                 ylabel('y');
%                 title('PumpL');
                
            else
                refreshdata(mesh_psiL,'caller');
                refreshdata(mesh_psiR,'caller');
                refreshdata(mesh_psi_LR,'caller');
                refreshdata(mesh_phi_L,'caller');
                refreshdata(mesh_phi_R,'caller');
%                 refreshdata(mesh_pumpL,'caller');
                
                if Pulse_add_flag==1
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_L_plot,{st;'|\psi_L|^2'});
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
    
    %Incoherent pulse
    if Pulse_add_flag==0 %stage 1
        if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
            fprintf('               pulse begins\n');
            Pulse_add_flag=1;%enter state 2
        end
    elseif Pulse_add_flag==1 %stage 2
            Pmp.time=time;
            Pse.time=double(time-Pse.force_pulse_time);
            adiabatic_term=0.5*(1+tanh(Pse.ada_coe*(Pse.time-Pse.ada_shift)));
            PumpL=abs(sqrt(Pmp.PL)*Pmp.u+sqrt(adiabatic_term)*sqrt(Pse.PL)*Pse.u).^2;
            PumpR=abs(sqrt(1-Pmp.PL)*Pmp.u+sqrt(adiabatic_term)*sqrt(1-Pse.PL)*Pse.u).^2;
            if M.GPU_calculation==2
               PumpL_gpu=gpuArray(PumpL);
               PumpR_gpu=gpuArray(PumpR); 
            end
        if (Pse.force_pulse_time+Pse.pulse_totaltime)<=time %stage 3
            fprintf('               pulse ends\n')
            Pmp.time=time;
            Pmp.Update_u();
            PumpL=Pmp.PL*abs(Pmp.u).^2;
            PumpR=(1-Pmp.PL)*abs(Pmp.u).^2;
            if M.GPU_calculation==2
               PumpL_gpu=gpuArray(PumpL);
               PumpR_gpu=gpuArray(PumpR); 
            end
            Pmp.clear_all();
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
if numlabs==1 %no SPMD
fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','zbuffer','Visible','off');

subplot(2,3,1);
mesh(x,y,abs(psiL).^2)
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi_L|^2');

subplot(2,3,2)
mesh(x,y,abs(psiR).^2);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title({file_name,'|\psi_R|^2'});

subplot(2,3,3)
mesh(x,y,abs(psiL).^2+abs(psiR).^2);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('|\psi_L|^2+|\psi_R|^2');

subplot(2,3,4)
mesh(x,y,wrapToPi(angle(psiL)));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi_L');

subplot(2,3,5)
mesh(x,y,wrapToPi(angle(psiR)));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi_R');

subplot(2,3,6)
plot(Psi.Time,Psi.Energy);
axis square;
xlabel('t');
title('Energy');

saveas(fig_E,strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'),'png');

else %for SPMD
    fig_E=figure('position',[100 100 500 500],'PaperPositionMode','auto','renderer','painter','Visible','off');
    plot(Psi.Time,Psi.Energy,'linewidth',2);
    set(gca,'fontsize',Psi.font_size);
    axis square;
    xlabel('t');
    title('Energy');
    print(fig_E,'-painter','-dpng','-r200',strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'));
end

close(fig_E);

fprintf('Storing Data\n ');

Time=Psi.Time;
Energy=Psi.Energy;
mod_psiLR_square=Psi.mod_psiLR_square;
if Psi.record_sxyz_int~=0
    sx_int=Psi.sx_int;
    sy_int=Psi.sy_int;
    sz_int=Psi.sz_int;
end
FS_fixed_int_r=Psi.fixed_int_r;
FS_int_tolerant=Psi.int_tolerant;
FS_int_r=Psi.int_r;

Pmp.clear_all();
Pse.clear_all();
Psi.clear_all();
M.clear_all();

Ori_Pump_Parameter=Pmp;
if Psi.record_finial_state==0 && Psi.record_psi~=0
    record_psi_t=Psi.record_psi_t;
    psiL_trim=Psi.psiL_trim;
    psiR_trim=Psi.psiR_trim;
    save(strcat(currentFolder,'\',dataFolder,'\',file_name,'.mat'),'FS_fixed_int_r','FS_int_tolerant','FS_int_r','Ori_Pump_Parameter','record_psi_t','psiL_trim','psiR_trim','Time','Energy','-mat','-v7.3');
elseif Psi.record_finial_state~=0
    Finial_State=complex(0,zeros(N,N,4));
    Finial_State(:,:,1)=psiL;
    Finial_State(:,:,2)=psiR;
    Finial_State(:,:,3)=nL;
    Finial_State(:,:,4)=nR;
    if Psi.record_sxyz_int~=0
        save(strcat(currentFolder,'\',dataFolder,'\','FinalState_',file_name,'.mat'),'FS_fixed_int_r','FS_int_tolerant','FS_int_r','Ori_Pump_Parameter','mod_psiLR_square','Time','Energy','sx_int','sy_int','sz_int','Finial_State','-mat','-v7.3');
    else
        save(strcat(currentFolder,'\',dataFolder,'\','FinalState_',file_name,'.mat'),'Ori_Pump_Parameter','mod_psiLR_square','Time','Energy','Finial_State','-mat','-v7.3');
    end
end

fprintf('\nSimulation function finished.\n\n')
toc;

end



