function     Toroidal_O_Single_GaussianPulse_FollowUp_Inco_Function(Pmp,Pse,Psi,M)

%Full simulation program with Incoherent pulse (Object oriented)
%save Time Energy
%record mov
%record trimmed psi
%Follow program with Incoherent pulse (Object oriented)
%save Time Energy
%record mov
%record trimmed psi
fprintf('\nSingle component followUp program.\n')
tic;

%set dataFolder name
currentFolder = pwd;
folder_name_parameter=sprintf('XYmax=%d_J=%.1f',Pmp.XYmax,Psi.J);
dataFolder=strcat('TorPsi_Single_',M.Simulation_version,'_IncoPulse_',folder_name_parameter);

load(strcat(currentFolder,'\',dataFolder,'\','Tor_Standard_Steady_State_Single_',...
    M.Simulation_version,'_m0_',num2str(Psi.Ini_m),'_l',num2str(Pmp.l),'_p',num2str(Pmp.p),'.mat'));
%Check parameters insistent
if ~Pmp.Check_Parameter(Ori_Pump_Parameter);
    Ori_Pump_Parameter
    fprintf('\nError: Parameters do not fit.\n\n');
    return;
end

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
J=Pmp.J; %general coefficient

%initialize psi
Psi.initialize();
Psi.step_E=length(Time);
Psi.Time(1,1:Psi.step_E)=double(Time);
Psi.Energy(1,1:Psi.step_E)=Energy;
Psi.step_E=Psi.step_E+1;
psi=Finial_State(:,:,1);
n=Finial_State(:,:,2);
step_0=floor(max(Time)/dt);

%construct Pump
Pmp.Build_u();
Pump=abs(Pmp.u).^2;

%construct Pulse
Pse.Build_u();
Pulse_add_flag=0; %0:pump only, 1:pump+pulse, 2:pulse ended

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

%Saved data file name
file_name_parameter=sprintf('J=%.1f_m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_T=%d',J,Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,floor(TotalTime));
Mov_file_name_psi=strcat('TorPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
Mov_file_name_sxyz=strcat('Torsxyz_',M.Simulation_version,'_IncoPulse_',file_name_parameter);
file_name=strcat('TorPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter);

if Pse.force_pulse_time<=TotalTime
    PT=min(floor(Pse.pulse_totaltime),floor(abs(TotalTime-Pse.force_pulse_time)));%Pulse Total Time
    file_name_parameter=sprintf('J=%.1f_m0=%d_Pmp0=%.1f_l=%d_p=%d_w0=%d_Pse0=%.1f_x0=%d_y0=%d_Sgm_x=%d_Sgm_y=%d_OmeP=%.2f_PT=%d_T=%d',...
                                 J,Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,Pse.Pbar0,Pse.x0,Pse.y0,Pse.Sigma_x,Pse.Sigma_y,Pse.Omega,PT,floor(TotalTime));
    Mov_file_name_psi=strcat('TorPsi_',M.Simulation_version,'_IncoPulse_',file_name_parameter,'_LG');
    Mov_file_name_sxyz=strcat('Torsxyz_',M.Simulation_version,'_IncoPulse_',file_name_parameter,'_LG');
    file_name=strcat('Psi_',M.Simulation_version,'_IcoPse_',file_name_parameter,'_LG');
end

%GPU calculation initializing
if M.GPU_calculation~=0
    GPU_Device = gpuDevice;
    gpuDevice(1);
    kinetic_factor_quarter_gpu=gpuArray(kinetic_factor_quarter);
    damping_gpu=gpuArray(damping);
    if M.GPU_calculation==2 %full GPU calculation
        psi_gpu=gpuArray(psi);
        n_gpu=gpuArray(n);
        Pump_gpu=gpuArray(Pump);
    end
    fprintf('\nGPU calculation enabled\n');
end

if floor(step_0*dt)>Pse.force_pulse_time
    fprintf('\nError: Force pulse time is smaller then the steady state time\n\n');
    return;
end

fftw('planner','exhaustive');

fprintf('\nCalculation begins\n\n');

for step_i=step_0:floor(TotalTime/dt)
    
    time=double(step_i*dt);
    
    %Incoherent pulse
    if Pulse_add_flag==0 %stage 1
        if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
            fprintf('               pulse begins\n');
            Pulse_add_flag=1;%enter state 2
        end
    end
    if Pulse_add_flag==1 %stage 2
        if Pmp.Omega~=0
            Pmp.time=time;
            Pmp.Update_u();
        end
        Pse.time=double(time-Pse.force_pulse_time);
        adiabatic_term=0.5*(1+tanh(Pse.ada_coe*(Pse.time-Pse.ada_shift)));
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
    
    %Core calculation----------------------------------
    if M.GPU_calculation==2 %Full GPU calculation
        
        n_gpu=n_gpu.*exp(-(GammaR+R*abs(psi_gpu).^2)*dt)+Pump_gpu*dt;
        psi1_gpu=ifft2((kinetic_factor_quarter_gpu.*fft2(psi_gpu)));
        psi1_gpu=psi1_gpu.*(exp((ua*abs(psi1_gpu).^2+gR*n_gpu+0.5i*(R*n_gpu-GammaC))*(-1i)*dt));
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
    
    if Pulse_add_flag==1
        fprintf('t = %.3f      Pulse added\n',time);
    else   fprintf('t = %.3f\n',time);
    end
    
    if M.GPU_calculation==2
        if (mod(step_i,Psi.E_step_interval)==0)...
                || (Psi.record_psi~=0 && time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0))...
                || (M.record_movie~=0 && time>=M.re_mov_begin && time<=M.re_mov_end && (mod(step_i,floor(M.re_mov_dt/dt))==0))
            psi=gather(psi_gpu);
        end
    end
    
    %record Time Energy
    if mod(step_i,Psi.E_step_interval)==0
        Psi.step_i=step_i;
        Psi.psi=psi;
        Psi.record_Time_Energy();
    end
    
    %store trimmed psi
    if Psi.record_psi~=0
        if time>=Psi.re_psi_begin && time<=Psi.re_psi_end && (mod(step_i,floor(Psi.re_psi_dt/dt))==0)
            Psi.psi=psi;
            Psi.record_Psi();
        end
    end
    
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
                
                pump_plot=subplot(2,2,2);
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

%clean up movie writer
if M.record_movie~=0
    if ~isempty(M.vidObjpsi)
        close(M.vidObjpsi);
        close(psi_fig);
    end
end
%clean up GPU calculation parameters
if M.GPU_calculation~=0
    if M.GPU_calculation==2
        psi=gather(psi_gpu);
        n=gather(n_gpu);
    end
    delete(GPU_Device);
end

%draw final state and Energy evolution
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

subplot(2,2,2)
mesh(x,y,Pump);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('Pump');

subplot(2,2,4)
plot(Psi.Time,Psi.Energy);
axis square;
title('Energy');

saveas(fig_E,strcat(currentFolder,'\',dataFolder,'\',file_name,'.png'),'png');
close(fig_E);

fprintf('Storing Data\n ');

Time=Psi.Time;
Energy=Psi.Energy;
Pmp.clear_all();
Pse.clear_all();
Psi.clear_all();
M.clear_all();

Ori_Pump_Parameter=Pmp;
if  Psi.record_psi~=0
    record_psi_t=Psi.record_psi_t;
    psi_trim=Psi.psi_trim;
    save(strcat(currentFolder,'\',dataFolder,'\',file_name,'.mat'),'Ori_Pump_Parameter','record_psi_t','psi_trim','Time','Energy','-mat','-v7.3');
else
    Finial_State=complex(0,zeros(N,N,2));
    Finial_State(:,:,1)=psi;
    Finial_State(:,:,2)=n;
    save(strcat(currentFolder,'\',dataFolder,'\','FinalState_',file_name,'.mat'),'Ori_Pump_Parameter','Time','Energy','Finial_State','-mat','-v7.3');
end

fprintf('\nSimulation function finished.\n\n')
toc;

end



