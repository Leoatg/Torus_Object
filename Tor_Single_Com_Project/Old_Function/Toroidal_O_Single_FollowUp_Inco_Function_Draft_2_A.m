function     Toroidal_O_Single_FollowUp_Inco_Function_Draft_2_A(Pmp,Pse,Psi,M)

%FollowUp simulation program with Incoherent pulse (Object oriented)
%record Time & Energy within Psi; or within GPU (Full GPU Simulation mode)
%record mov
%record trimmed psi

fprintf('\nSingle component followUp program.\n')
tic;

%% load data
%set dataFolder name
currentFolder = pwd;
folder_name_parameter=sprintf('XYmax=%d',Pmp.XYmax);
dataFolder=strcat('TorPsi_Single_',M.Simulation_version,'_IncoPulse_',folder_name_parameter);

% load(strcat(currentFolder,'\',dataFolder,'\','Tor_Standard_Steady_State_Single_',...
%     M.Simulation_version,'_m0_',num2str(Psi.Ini_m),'_l',num2str(Pmp.l),'_p',num2str(Pmp.p),'.mat'));

load(strcat(currentFolder,'\',dataFolder,'\','FinalState_Tor_Pert_GmaRr=1.4_m0=66_Pmp0=2.5_l=5_p=0_w0=28_PerS=0.1_PerF=-46_PerT=80_XYmax=90_N=1024_T=1200.mat'));

%Check parameters insistent
if ~Pmp.Check_Parameter(Ori_Pump_Parameter);
    Ori_Pump_Parameter
    fprintf('\nError: Parameters do not fit.\n\n');
    return;
end

%% initialize parameters
N=Pmp.N;
XYmax=Pmp.XYmax;

TotalTime=Pmp.TotalTime;
dt=Pmp.dt;
if TotalTime<=max(Time);
    fprintf('\nError: TotalTime is smaller then steady state time\n\n');
    return;
end

%Coefficients in the ODGPE: (parameters from the PRA)
ua=Pmp.ua; %same spin polariton-polariton scattering
gR=Pmp.gR; %LP interaction with incoherent reservoir
GammaC=Pmp.GammaC;
GammaR=Pmp.GammaR; % GammaR_ratio=GammaR/GammaC
R=Pmp.R; %stimulated scattering rate

%construct Pump
Pmp.Build_u();
Pump=abs(Pmp.u).^2;

%construct Pulse
Pse.Build_u();
Pulse_add_flag=0; %0:pump only, 1:pump+pulse, 2:pulse ended
pulse_Omega=Pse.Omega;
Pse.Omega=0;
pulse_totate_begin_t=0;

%construct lattice size and k vector space
x=Pmp.x;
y=Pmp.y;
kvector=fftshift(-N/2:N/2-1)*2*pi/(2*XYmax);
[kx, ky]=meshgrid(kvector, kvector);
kinetic_factor_quarter=exp((-1i*(kx.^2+ky.^2)*dt)/4);

%Absorbing boundaries
r_damping=XYmax-Psi.damping_shift;%XYmax-14; %radius for top-hat damping
a_damping=Psi.a_damping; %slope for the damping mesa
damping=0.25*(1+tanh(a_damping*((x.^2+y.^2).^0.5+r_damping))).*(1+tanh(a_damping*(-(x.^2+y.^2).^0.5+r_damping)));
Psi.int_r=FS_fixed_int_r;

%initialize psi
Psi.initialize();
Psi.step_E=length(Time);
Psi.Time(1,1:Psi.step_E)=double(Time);
Psi.Energy(1,1:Psi.step_E)=Energy;
Psi.Lz(1,1:Psi.step_E)=Lz;
Psi.step_E=Psi.step_E+1;
psi=Finial_State(:,:,1);
n=Finial_State(:,:,2);
step_0=floor(max(Time)/dt);

%% set data file name
%Saved data file name
file_name_parameter=sprintf('m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_T=%d',Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,floor(TotalTime));
Mov_file_name_psi=strcat('TorPsi_IncoPulse_',file_name_parameter);
Mov_file_name_sxyz=strcat('Torsxyz_IncoPulse_',file_name_parameter);
file_name=strcat('TorPsi_IncoPulse_',file_name_parameter);

if Pse.force_pulse_time<=TotalTime
    %PT=min(floor(Pse.pulse_totaltime),floor(abs(TotalTime-Pse.force_pulse_time)));%Pulse Total Time
    file_name_parameter=sprintf('m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_Pse0=%.2f_LP=%d_pP=%d_w0P=%d_OmeP=%.2f_PWT=%d_PseT=%d_T=%d',Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,Pse.Pbar0,Pse.l,Pse.p,Pse.w0,pulse_Omega,Pse.pulse_rotate_t,floor(Pse.pulse_totaltime),floor(TotalTime));
    Mov_file_name_psi=strcat('TorPsi_IncoPulse_',file_name_parameter,'_LG');
    Mov_file_name_sxyz=strcat('Torsxyz_IncoPulse_',file_name_parameter,'_LG');
    file_name=strcat('Psi_IcoPse_',file_name_parameter,'_LG');
end

%% GPU calculation initializing
if M.GPU_calculation~=0
    GPU_Device = gpuDevice;
    gpuDevice(1);
    kinetic_factor_quarter_gpu=gpuArray(kinetic_factor_quarter);
    damping_gpu=gpuArray(damping);
    if M.GPU_calculation==2 %full GPU calculation
        psi_gpu=gpuArray(psi);
        n_gpu=gpuArray(n);
        Pump_gpu=gpuArray(Pump);
        dt_gpu=gpuArray(dt);
        ua=gpuArray(ua);
        gR=gpuArray(gR);
        R=gpuArray(R);
        hspace=gpuArray(Psi.hspace);
        GammaC=gpuArray(GammaC);
        GammaR=gpuArray(GammaR);
        Energy_gpu=gpuArray(Psi.Energy);
        Lz_gpu=gpuArray(Psi.Lz);
        step_E_gpu=gpuArray(Psi.step_E);
        int_area_gpu=gpuArray(Psi.int_area);
        x_gpu=gpuArray(x);
        y_gpu=gpuArray(y);
        step_E=Psi.step_E;
        clear Energy Lz
    end
    fprintf('\nGPU calculation enabled\n');
end

if floor(step_0*dt)>Pse.force_pulse_time
    fprintf('\nError: Force pulse time is smaller then the steady state time\n\n');
    return;
end

fftw('planner','exhaustive');

fprintf('\nCalculation begins\n\n');

%% mian calculation
for step_i=step_0:floor(TotalTime/dt)
    
    time=double(step_i*dt);
    
    %% update pulse
    %Incoherent pulse
    if Pulse_add_flag==0 %stage 1
        if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
            fprintf('               pulse begins\n');
            Pulse_add_flag=1;%enter state 2
        end
    end
    if Pulse_add_flag==1 %stage 2
        if Pmp.Omega~=0 %update pump
            Pmp.time=time;
            Pmp.Update_u();
        end
        if (time-Pse.force_pulse_time)>=Pse.pulse_rotate_t
            Pse.Omega=pulse_Omega;
            if pulse_totate_begin_t==0;
            pulse_totate_begin_t=double(time-Pse.force_pulse_time);
            end
        end
        pulse_adiabatic_t=double(time-Pse.force_pulse_time);
        Pse.time=pulse_adiabatic_t-pulse_totate_begin_t;
        adiabatic_term=0.5*(1+tanh(Pse.ada_coe*(pulse_adiabatic_t-Pse.ada_shift)));
        Pse.Update_u();
        Pump=abs(Pmp.u+sqrt(adiabatic_term)*Pse.u).^2;
        if M.GPU_calculation==2
            Pump_gpu=gpuArray(Pump);
        end
        if (Pse.force_pulse_time+Pse.pulse_totaltime)<=time %stage 3
            fprintf('               pulse ends\n')
            Pmp.time=time;
            Pmp.Update_u();
            Pump=abs(Pmp.u).^2;
            if M.GPU_calculation==2
               Pump_gpu=gpuArray(Pump);
            end
            Pmp.clear_all();
            Pse.clear_all();
            Pulse_add_flag=2;%enter state 3
        end
    end %pulse update end
    
    %%    
    %Core calculation----------------------------------
    if M.GPU_calculation==2 %Full GPU calculation
        
        n_gpu=n_gpu.*exp(-(GammaR+R*abs(psi_gpu).^2)*dt_gpu)+Pump_gpu*dt_gpu;
        psi1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psi_gpu)));
        psi1_gpu=psi1_gpu.*(exp((ua*abs(psi1_gpu).^2+gR*n_gpu+0.5i*(R*n_gpu-GammaC))*(-1i)*dt_gpu));
        psi_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psi1_gpu))).*damping_gpu;
        
    else %CPU or Half GPU
        
        %reservoirs update
        n=n.*exp(-(GammaR+R*abs(psi).^2)*dt)+Pump*dt;
        
        %advance spacially the first time
        if M.GPU_calculation==1 %GPU FFT
            psi1=gather(ifft2(kinetic_factor_quarter_gpu.*fft2(gpuArray(psi))));
        else
            psi1=ifft2((kinetic_factor_quarter.*fft2(psi)));
        end
        
        %update psiL1,R1 by the pseudo-potential
        psi1=psi1.*(exp((ua*abs(psi1).^2+gR*n+0.5i*(R*n-GammaC))*(-1i)*dt));
        
        %advance the second time with boundary condition (damping)
        if M.GPU_calculation==1 %GPU FFT
            psi=gather(ifft2((kinetic_factor_quarter_gpu.*fft2(gpuArray(psi1)))).*damping_gpu);
        else
            psi=ifft2((kinetic_factor_quarter.*fft2(psi1))).*damping;
        end
    end
    %----------------------------------
    
    %% status & data recording
    if Pulse_add_flag==1
        fprintf('t = %.3f      Pulse added\n',time);
    else   fprintf('t = %.3f\n',time);
    end
    
    Psi_data_updated_flag=0;
    
    if M.GPU_calculation==2
        if  (Psi.record_psi~=0 && time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0))...
                || (M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0))
            psi=gather(psi_gpu);
            n=gather(n_gpu);
        end
    end
    
    %record Time Energy
    if mod(step_i,Psi.E_step_interval)==0
        
        if M.GPU_calculation==2 %record Energy within the GPU
            
            [psix, psiy]=gradient(psi_gpu,hspace);
            mod_psi_square=abs(psi_gpu).^2;
            psi_int_gpu=hspace^2*sum(sum(mod_psi_square.*int_area_gpu));
            
            Lz_gpu(step_E_gpu)=hspace^2*sum(sum((-1i)*conj(psi_gpu).*(x_gpu.*psiy-y_gpu.*psix).*int_area_gpu))/psi_int_gpu;
 
            Energy_gpu(step_E_gpu)=(hspace^2*sum(sum((0.5*(abs(psix).^2+abs(psiy).^2)+gR*n_gpu.*mod_psi_square+0.5*ua*mod_psi_square.^2).*int_area_gpu)))/psi_int_gpu;
            
            step_E_gpu=step_E_gpu+1;
            
            Time(step_E)=double(dt*step_i);
            step_E=step_E+1;
            
        else %record Energy
        Psi.step_i=step_i;
        Psi.n=n;
        Psi.psi=psi;
        Psi_data_updated_flag=1;
        Psi.record_Time_Energy();%Lz included
        end
    end
    
    %store trimmed psi
    if Psi.record_psi~=0
        if time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0)
            if Psi_data_updated_flag==0
                Psi.psi=psi;
                Psi.n=n;
                Psi_data_updated_flag=1;
            end
            Psi.record_Psi();
        end
    end
    
    %% record Mov
    %record Psi Mov
    if M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0)
        
        mod_psi=abs(psi);
        mod_psi_square=mod_psi.^2;
        phi=wrapToPi(angle(psi));
        
        if M.record_mov_psi~=0
            
            if M.re_psi_initialized==0 %psi windows and mesh initialization
                
                psi_fig=figure('position',[100 200 900 700],'renderer','zbuffer');
                M.re_psi_initialized=1;
                M.mov_file_name_psi=Mov_file_name_psi;
                M.currentFolder=currentFolder;
                M.dataFolder=dataFolder;
                
                psi_plot=subplot(2,2,1);
                mesh_psi=mesh(x,y,mod_psi_square);
                axis square;
                set(mesh_psi,'ZDataSource','mod_psi_square');
                set(psi_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                if Pulse_add_flag==1%Pse.force_pulse_time>=time && (Pse.force_pulse_time+Pse.pulse_totaltime)<=time
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_plot,{st;'|\psi|^2'});
                
                psi_plot2=subplot(2,2,2);
                mesh_psi2=mesh(x,y,mod_psi_square);
                axis square;
                set(mesh_psi2,'ZDataSource','mod_psi_square');
                set(psi_plot2,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                zlabel('|\psi|^2');
                colorbar;
                
                phi_plot=subplot(2,2,3);
                mesh_phi=mesh(x,y,phi);
                axis square;
                set(mesh_phi,'ZDataSource','phi');
                set(phi_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('\phi');
                
                pump_plot=subplot(2,2,4);
                mesh_pump=mesh(x,y,Pump);
                axis square;
                set(mesh_pump,'ZDataSource','Pump');
                set(pump_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('Pump');
                
            else
                refreshdata(mesh_psi,'caller');
                refreshdata(mesh_psi2,'caller');
                refreshdata(mesh_phi,'caller');
                refreshdata(mesh_pump,'caller');
                
                if Pulse_add_flag==1
                    st=sprintf('Pulse added        t=%.1f',time);
                else   st=sprintf('                   t=%.1f',time);
                end
                title(psi_plot,{st;'|\psi|^2'});
            end
            M.Mov_psi(psi_fig);
        end % psi movie end
        
    end %record two movies end
    
end %end of main loop

%% clean up movie writer
if M.record_movie~=0
    if ~isempty(M.vidObjpsi)
        close(M.vidObjpsi);
        close(psi_fig);
    end
end

%% prepare for finial record
if M.GPU_calculation==2
    Energy=gather(Energy_gpu);
    Lz=gather(Lz_gpu);
    psi=gather(psi_gpu);
    n=gather(n_gpu);
    %Time is recored during Energy_gpu
else
    Time=Psi.Time;
    Energy=Psi.Energy;
    Lz=Psi.Lz;
end

%clean up GPU calculation parameters
if M.GPU_calculation~=0
    delete(GPU_Device);
end

%% draw final state and Energy evolution
if numlabs==1 %no SPMD
fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','zbuffer','Visible','off');

subplot(2,2,1);
mesh(x,y,abs(psi).^2);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi|^2');

subplot(2,2,3);
mesh(x,y,wrapToPi(angle(psi)));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi');

subplot(2,2,2);
plot(Time,real(Lz));
axis square;
set(gca,'fontsize',Psi.font_size);
xlabel('t');
ylabel('L_z');
st=sprintf('                                   L_z final = %.1f',real(Lz(length(Lz))));
title(st);

subplot(2,2,4);
plot(Time,Energy);
axis square;
set(gca,'fontsize',Psi.font_size);
xlabel('t');
title('Energy (normalized)');

saveas(fig_E,strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'),'png');
else %for SPMD
    fig_E=figure('position',[100 100 800 800],'PaperPositionMode','auto','renderer','painter','Visible','off', 'PaperType', 'a4');
    
    plot(Time,Psi.Energy_unmal,'linewidth',2);
    set(gca,'fontsize',Psi.font_size);
    axis square;
    xlabel('t');
    title('Energy (normalized)');
    
    print(fig_E,'-dpng',strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'));
end

close(fig_E);

fprintf('Storing Data\n ');

FS_fixed_int_r=Psi.int_r;
Pmp.clear_all();
Pse.clear_all();
Psi.clear_all();
M.clear_all();

%% save data to mat file
Ori_Pump_Parameter=Pmp;
Finial_State=complex(0,zeros(N,N,2));
Finial_State(:,:,1)=psi;
Finial_State(:,:,2)=n;
if  Psi.record_psi~=0 %time series of psi is saved
    record_psi_t=Psi.record_psi_t;
    psi_trim=Psi.psi_trim;
    n_trim=Psi.n_trim;
    save(strcat(currentFolder,'\',dataFolder,'\',file_name,'.mat'),'Ori_Pump_Parameter','Finial_State','FS_fixed_int_r','record_psi_t','psi_trim','n_trim','Time','Energy','Lz','-mat','-v7.3');
else                 %save only the finial state
    save(strcat(currentFolder,'\',dataFolder,'\','FinalState_',file_name,'.mat'),'Ori_Pump_Parameter','Finial_State','FS_fixed_int_r','Time','Energy','Lz','-mat','-v7.3');
end

fprintf('\nSimulation function finished.\n\n')
toc;

end



