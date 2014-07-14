function     Toroidal_O_Single_Pert_Simulation_Inco_Function_try2(Pmp,Pse,Psi,M)

%Full simulation program with Incoherent pulse (Object oriented)
%save Time Energy
%record mov
%record trimmed psi
fprintf('\nSingle component Pert full simulation program.\n')

N=Pmp.N;
XYmax=Pmp.XYmax;
global OS;
global NoDisplay_Flag;

TotalTime=Pmp.TotalTime;
dt=Pmp.dt;

%% Coefficients in the ODGPE: (parameters from the PRA)
ua=Pmp.ua; %same spin polariton-polariton scattering
gR=Pmp.gR; %LP interaction with incoherent reservoir
GammaC=Pmp.GammaC;
GammaR=Pmp.GammaR; % GammaR_ratio=GammaR/GammaC
R=Pmp.R; %stimulated scattering rate
hspace=Pmp.hspace;

%% Construct Pump
Pmp.Build_u();
Pump=abs(Pmp.u).^2; % Intensity
R0=Pmp.w0*sqrt(Pmp.l/2);

%calculate effective Pump power P_eff
N_mid=floor(N/2+1);
P_th=Pmp.P0;
Pth_cutoff=1.45;
P_eff=-1;
for it=N_mid:N
    if Pump(it,N_mid)>=(P_th*Pth_cutoff)
        N1_P_eff=it;
        break;
    end
end
for it=(N1_P_eff+1):N
    if Pump(it,N_mid)<=((P_th*Pth_cutoff))
        N2_P_eff=it;
        break;
    end
end
if N2_P_eff>N1_P_eff
    r_b=(N1_P_eff-N_mid)*hspace;
    r_e=(N2_P_eff-N_mid)*hspace;
    P_eff=(sum(Pump(N1_P_eff:N2_P_eff,N_mid))/length(N1_P_eff:N2_P_eff))/P_th;
end

%% Construct Pulse
Pse.Build_u();
Pulse_add_flag=0; %stage: 0:pump only, 1:pump+pulse, 2:pulse ended

%% Construct lattice size and k vector space
x=Pmp.x;
y=Pmp.y;
kvector=fftshift(-N/2:N/2-1)*2*pi/(2*XYmax);
[kx, ky]=meshgrid(kvector, kvector);
kinetic_factor_quarter=exp((-1i*(kx.^2+ky.^2)*dt)/4);

%% Absorbing boundaries
r_damping=XYmax-Psi.damping_shift;%radius for top-hat damping
a_damping=Psi.a_damping; %slope for the damping mesa
damping=0.25*(1+tanh(a_damping*((x.^2+y.^2).^0.5+r_damping))).*(1+tanh(a_damping*(-(x.^2+y.^2).^0.5+r_damping)));
%x_axis=-XYmax:Pmp.hspace:XYmax;
%set int_r
for it=floor(N/2+1):N
    if damping(it,floor(N/2+1))<=0.95
       Psi.int_r=(it-floor(N/2+1))*Pmp.hspace;
       break;
    end
end

%% Initialize psi
Psi.initialize();
psi=Psi.psi;
n=Psi.n; %reservoir

%% Saved data file name
file_name_parameter=sprintf('GmaRr=%.1f_m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_PerS=%g_PerF=%.1f_PerT=%d_XYmax=%d_N=%d_T=%d',Psi.GammaR/Psi.GammaC,Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,M.Pert_Strength,M.Pert_Fre,M.Pert_time,Pmp.XYmax,N,floor(TotalTime));
Mov_file_name_psi=strcat('Tor_Pert_',file_name_parameter);
file_name=strcat('Tor_Pert_',file_name_parameter);

if Pse.force_pulse_time<=TotalTime
    file_name_parameter=sprintf('m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_Pse0=%.1f_LP=%d_pP=%d_w0P=%d_OmegaP=%.4f_N=%d_T=%d',Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,Pse.Pbar0,Pse.l,Pse.p,Pse.w0,Pse.Omega,N,floor(TotalTime));
    Mov_file_name_psi=strcat('Tor_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
    file_name=strcat('Psi_Pert_',file_name_parameter);
end

currentFolder = pwd;
folder_name_parameter=sprintf('XYmax=%d',Pmp.XYmax);
dataFolder=strcat('TorPsi_Single_',M.Simulation_version,'_IncoPulse_',folder_name_parameter);

if exist(dataFolder,'dir')==0
    mkdir(dataFolder);
end

%% GPU calculation initializing
if M.GPU_calculation~=0
    
    fprintf('\nGPU calculation enabled\n\n');
    if M.GPU_calculation>=3
        GPU_Device = gpuDevice(M.GPU_calculation-2);
    else
        GPU_Device = gpuDevice;
    end
    %gpuDevice(1);
    fprintf(['GPU Device: ',GPU_Device.Name,' selected.\n']);
    kinetic_factor_quarter_gpu=gpuArray(kinetic_factor_quarter);
    damping_gpu=gpuArray(damping);
    if M.GPU_calculation>=2 %full GPU calculation
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
        
        Energy_unmal_gpu=Energy_gpu;
        psi_int_gpu=Energy_gpu;
    end
end

%% FFT prepare
fftw('planner','exhaustive');

fprintf('\nCalculation begins\n\n');

%% main calculation
for step_i=1:floor(TotalTime/dt)
    
    time=double(step_i*dt);
    
    %cool down Tesla
    if M.GPU_calculation==3 && ((TotalTime-time)>1)
        if N>=2048
            if mod(time,2)==0 && (time>1)
                fprintf('               Tesla cool down period\n');
                if time>20
                    pause(16);
                else
                    pause(12);
                end
            end
        elseif N>1024
            if mod(time,5)==0 && (time>6)
                fprintf('               Tesla cool down period\n');
                pause(20);
            end
        end
    end
    
    %% Initialize Perturbation
    if M.Pert_time>0 && M.Pert_Strength~=0 && time>M.Pert_time  %nagetive value means perturbation added
        if M.GPU_calculation>=2
            psi=gather(psi_gpu);
        end
        psi_mode=abs(psi);     
        phi0=Psi.Ini_m; % steady state m
        theta=atan2(y,x);
        k_per=M.Pert_Fre;
        rho_1=-1; %u_0
        rho_2=1; %v_0
        
        Pert=rho_1*exp(1i*(phi0+k_per)*theta)+rho_2*exp(1i*(phi0-k_per)*theta);
        %Pert=Pert+exp(1i*(phi0+1.5*k_per+pi/3)*theta)+exp(1i*(phi0-1.5*k_per-pi/3)*theta);
        %Pert=Pert+0*(rand(N)-1);
        Pert=psi_mode.*Pert*M.Pert_Strength;

        if M.GPU_calculation>=2
           psi_gpu=gpuArray(psi+Pert); 
        else
           psi=psi+Pert;
        end
        fprintf(sprintf('               perturbation added at t=%.1f\n',time));
        M.Pert_time=-1;
        
    end
    
    %% update pulse    
    %Incoherent pulse
    if Pulse_add_flag==0 %stage 1
        if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
            fprintf('               pulse begins\n');
            Pulse_add_flag=1;%enter state 2
        end
    end
    if Pulse_add_flag==1 %stage 2
        Pmp.time=time;
        if Pmp.Omega~=0
            Pmp.Update_u();
        end
        Pse.time=time-Pse.force_pulse_time;
        Pse.Update_u();
        Pump=abs(Pmp.u+Pse.u).^2;
        if M.GPU_calculation>=2
            Pump_gpu=gpuArray(Pump);
        end
        if (Pse.force_pulse_time+Pse.pulse_totaltime)<=time %stage 3
            fprintf('               pulse ends\n')
            Pmp.time=time;
            Pmp.Update_u();
            Pump=abs(Pmp.u).^2;
            if M.GPU_calculation>=2
                Pump_gpu=gpuArray(Pump);
            end
            Pmp.clear_all();
            Pse.clear_all();
            Pulse_add_flag=2;%enter state 3
        end
    end %pulse update end
    
    %%     
    %Core calculation----------------------------------
    if M.GPU_calculation>=2 %Full GPU calculation
        
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
    if OS==0
        if Pulse_add_flag==1
            fprintf('t = %.4f      Pulse added\n',time);
        else   fprintf('t = %.4f\n',time);
        end
    end
    
    Psi_data_updated_flag=0;
    
    if M.GPU_calculation>=2
        if  (Psi.record_psi~=0 && time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0))...
                || (M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0))
            psi=gather(psi_gpu);
            n=gather(n_gpu);
        end
    end
    
    %record Time Energy
    if mod(step_i,Psi.E_step_interval)==0
        
        if M.GPU_calculation>=2 %record Energy within the GPU
            
            [psix, psiy]=gradient(psi_gpu,hspace);
            mod_psi_square=abs(psi_gpu).^2;
            psi_int_gpu(step_E_gpu)=hspace^2*sum(sum(mod_psi_square.*int_area_gpu));
            
            Lz_gpu(step_E_gpu)=hspace^2*sum(sum((-1i)*conj(psi_gpu).*(x_gpu.*psiy-y_gpu.*psix).*int_area_gpu))/psi_int_gpu(step_E_gpu);
            
            Energy_unmal_gpu(step_E_gpu)=(hspace^2*sum(sum((0.5*(abs(psix).^2+abs(psiy).^2)+gR*n_gpu.*mod_psi_square+0.5*ua*mod_psi_square.^2).*int_area_gpu)));
            
            Energy_gpu(step_E_gpu)=Energy_unmal_gpu(step_E_gpu)/psi_int_gpu(step_E_gpu);
            
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
            end
            Psi.record_Psi();
        end
    end
    
    %% record Mov    
    %record Psi Mov
    if (M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0) )|| ((M.record_movie~=0)&&((step_i==1)&& M.re_mov_begin==0))
        
        mod_psi=abs(psi);
        mod_psi_square=mod_psi.^2;
        phi=wrapToPi(angle(psi));
        
        if M.record_mov_psi~=0
            
            if M.re_psi_initialized==0 %psi windows and mesh initialization
                
                psi_fig=figure('position',[80 80 1000 700],'renderer','zbuffer');
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
                
                n_plot=subplot(2,2,4);
                mesh_n=mesh(x,y,n);
                axis square;
                set(mesh_n,'ZDataSource','n');
                set(n_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                xlabel('x');
                ylabel('y');
                view(0,90);
                zoom(Psi.zoom_factor);
                colorbar;
                title('n');
                
            else
                refreshdata(mesh_psi,'caller');
                refreshdata(mesh_psi2,'caller');
                refreshdata(mesh_phi,'caller');
                refreshdata(mesh_n,'caller');
                
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
if M.GPU_calculation>=2
    Energy=gather(Energy_gpu);
    Lz=gather(Lz_gpu);
    psi=gather(psi_gpu);
    n=gather(n_gpu);
    
    Energy_unmal=gather(Energy_unmal_gpu);
    Psi_int=gather(psi_int_gpu);
    %Time is recored during Energy_gpu
else
    Time=Psi.Time;
    Energy=Psi.Energy;
    Lz=Psi.Lz;
    
    Energy_unmal=Psi.Energy_unmal;
    Psi_int=Psi.Psi_int;
end

%clean up GPU calculation parameters
if M.GPU_calculation~=0
    delete(GPU_Device);
end

%% draw final state and Energy evolution
if numlabs==1 && OS==0 && NoDisplay_Flag==0 %no SPMD, windows, Display
fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','zbuffer','Visible','off');

subplot(2,2,1);
mesh(x,y,abs(psi).^2)
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi|^2');

subplot(2,2,3)
mesh(x,y,wrapToPi(angle(psi)));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi');

% subplot(2,2,2)
% mesh(x,y,Pump);
% pos_tem=get(gca,'position');
% pos_tem(1)=0.3692;
% pos_tem2=pos_tem;
% set(gca,'position',pos_tem);
% axis square;
% set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
% view(0,90);
% zoom(Psi.zoom_factor);
% colorbar;
% xlabel('x');
% %ylabel('y');
% title('Pump');


subplot(2,2,2)
plot(Time,Energy_unmal,'linewidth',2);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
pos_tem2=pos_tem;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size);
xlabel('t');
%ylabel('y');
title('Energy (Unnormalized)');

subplot(2,2,4);
mesh(x,y,n);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('n');

%right bottom corner
axes('position',[ 0.6420    0.1256    0.3347    0.3412]);
plot(Time,Energy,'linewidth',2);
pos_tem(1)=0.5859;
set(gca,'position',pos_tem);
set(gca,'fontsize',Psi.font_size);
if max(num2str(get(gca,'YTick')))<1
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'));
end
axis square;
xlabel('t');
title('Energy (normalized)');

%right upper corner
pos_tem2(1)=0.5859;
axes('position',pos_tem2);
plot(Time,real(Lz),'linewidth',2);
set(gca,'fontsize',Psi.font_size);
axis square;
xlabel('t');
title(sprintf('L_z (normalized)  m_0=%d;  L_z final=%.1f',Psi.Ini_m,real(Lz(length(Lz)))));

saveas(fig_E,strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'),'png');

else %for SPMD or Linux
    fprintf('Saving figures \n ');
    
%     fig_E=figure('position',[100 100 800 800],'PaperPositionMode','auto','renderer','painter','Visible','off');
%     plot(Psi.Time,Psi.Energy_unmal,'linewidth',2);
%     set(gca,'fontsize',Psi.font_size);
%     axis square;
%     xlabel('t');
%     title('Energy (un-normalized)');
    
fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','zbuffer','Visible','off');

subplot(2,2,1);
mesh(x,y,abs(psi).^2)
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
title('|\psi|^2');

subplot(2,2,3)
mesh(x,y,wrapToPi(angle(psi)));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi');

% subplot(2,2,2)
% mesh(x,y,Pump);
% pos_tem=get(gca,'position');
% pos_tem(1)=0.3692;
% pos_tem2=pos_tem;
% set(gca,'position',pos_tem);
% axis square;
% set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
% view(0,90);
% zoom(Psi.zoom_factor);
% colorbar;
% xlabel('x');
% %ylabel('y');
% title('Pump');


subplot(2,2,2)
plot(Time,Energy_unmal,'linewidth',2);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
pos_tem2=pos_tem;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size);
xlabel('t');
%ylabel('y');
title('Energy (Unnormalized)');

subplot(2,2,4);
mesh(x,y,n);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('n');

%right bottom corner
axes('position',[ 0.6420    0.1256    0.3347    0.3412]);
plot(Time,Energy,'linewidth',2);
pos_tem(1)=0.5859;
set(gca,'position',pos_tem);
set(gca,'fontsize',Psi.font_size);
if max(num2str(get(gca,'YTick')))<1
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'));
end
axis square;
xlabel('t');
title('Energy (normalized)');

%right upper corner
pos_tem2(1)=0.5859;
axes('position',pos_tem2);
plot(Time,real(Lz),'linewidth',2);
set(gca,'fontsize',Psi.font_size);
axis square;
xlabel('t');
title(sprintf('L_z (normalized)  m_0=%d;  L_z final=%.1f',Psi.Ini_m,real(Lz(length(Lz)))));



    print(fig_E,'-dppm','-r120',strcat(currentFolder,'\',dataFolder,'\',file_name,'.ppm'));

end
close(fig_E);

FS_fixed_int_r=Psi.int_r;
Pmp.clear_all();
Pse.clear_all();
Psi.clear_all();
M.clear_all();


%% save data to mat file
Ori_Pump_Parameter=Pmp;

if OS==0
    file_path=[currentFolder,'\',dataFolder,'\',file_name,'.mat'];
else
    file_path=[currentFolder,'/',dataFolder,'/',file_name,'.mat'];
end
if Psi.record_psi~=0 %time series of psi is saved
    record_psi_t=Psi.record_psi_t;
    psi_trim=Psi.psi_trim;
    n_trim=Psi.n_trim;
    save(file_path,'Finial_State','N2_P_eff','N1_P_eff','Pth_cutoff','P_eff','FS_fixed_int_r','Ori_Pump_Parameter','record_psi_t','psi_trim','n_trim','Time','Energy','Psi_int','Lz','-mat','-v7.3');
else
    if Psi.record_finial_state==0 %no data saved
        fprintf('No Data Saved.\n ');
    elseif Psi.record_finial_state==1 %save Pump and Energy
        fprintf('Storing Data\n ');
        save(file_path,'N2_P_eff','N1_P_eff','Pth_cutoff','P_eff','FS_fixed_int_r','Ori_Pump_Parameter','Time','Energy','Psi_int','Lz','-mat','-v7.3');
    elseif Psi.record_finial_state==2 %save the finial state
        fprintf('Storing Data\n ');
        Finial_State=complex(0,zeros(N,N,2));
        Finial_State(:,:,1)=psi;
        Finial_State(:,:,2)=n;
        save(file_path,'N2_P_eff','N1_P_eff','Pth_cutoff','P_eff','FS_fixed_int_r','Ori_Pump_Parameter','Time','Energy','Psi_int','Finial_State','Lz','-mat','-v7.3');
    end
end

fprintf('\nSimulation function finished.\n\n')

end



