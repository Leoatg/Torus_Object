
function Toroidal_TestPump_Function(Pmp,Pse,Psi,M)
%Test Pump and Pulse--LG modes

currentFolder = pwd;
dataFolder='test_pump';
if exist(dataFolder,'dir')==0
    mkdir(dataFolder);
end

%Saved data file name

file_name_parameter=sprintf('J=%.1f_Pmp0=%.2f_l=%d_p=%d_w0=%d_T=%d',Pmp.J,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,floor(Pmp.TotalTime));
Mov_file_name_psi=strcat('Test_Pump_',M.Simulation_version,'_',file_name_parameter);

if Pse.force_pulse_time<=Pmp.TotalTime
    if isa(Pse,'LG_pulse')
        file_name_parameter=sprintf('J=%.1f_m0=%d_Pmp0=%.2f_l=%d_p=%d_w0=%d_Pse0=%.2f_LP=%d_pP=%d_w0P=%d_OmegaP=%.5f_T=%d',Pmp.J,Psi.Ini_m,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,Pse.Pbar0,Pse.l,Pse.p,Pse.w0,Pse.Omega,floor(Pmp.TotalTime));
    elseif isa(Pse,'Gaussian_pulse')
        file_name_parameter=sprintf('J=%.1f_Pmp0=%.2f_l=%d_p=%d_w0=%d_Pse0=%.2f_OmegaP=%.5f_T=%d',Pmp.J,Pmp.Pbar0,Pmp.l,Pmp.p,Pmp.w0,Pse.Pbar0,Pse.Omega,floor(Pmp.TotalTime));
    end
    Mov_file_name_psi=strcat('Test_Pump_',M.Simulation_version,'_',file_name_parameter);
end

XYmax=Pmp.XYmax;

Pmp.Build_u();
Pse.Build_u();
Pulse_add_flag=0;

if Pmp.component_number==1
    Pump=abs(Pmp.u).^2;
    Pump_phase=wrapToPi(angle(Pmp.u));
    Pulse=abs(Pse.u).^2;
    Pulse_phase=wrapToPi(angle(Pse.u));
    Pump_combined=Pump;
    Pump_combined_phase=Pump_phase;
elseif Pmp.component_number==2
    
end

x=Pmp.x;
y=Pmp.y;

dt=Pmp.dt;
Step_Mov=floor(M.re_mov_dt/dt);

figure;
plot(Pump(:,Pmp.N/2));

Pump_and_Pulse_mesh=1;

for step_i=0:floor(Pmp.TotalTime/dt)
    
    if mod(step_i,Step_Mov)==0
        
        time=double(step_i*dt);
        
        %Incoherent pulse
        if Pulse_add_flag==0 %stage 1
            if Pse.force_pulse_time<=time && (Pse.force_pulse_time+Pse.pulse_totaltime)>=time
                fprintf('               pulse begins\n');
                Pulse_add_flag=1;%enter state 2
            end
        end
        if Pulse_add_flag==1 %stage 2
            Pmp.time=time;
            Pse.time=double(time-Pse.force_pulse_time);
            Pmp.Update_u();
            Pse.Update_u();
            adiabatic_term=1;%0.5*(1+tanh(Pse.ada_coe*(Pse.time-Pse.ada_shift)));
            if Pse.component_number==1
                Pump_combined=abs(Pmp.u+sqrt(adiabatic_term)*Pse.u).^2;
                Pump_combined_phase=wrapToPi(angle(Pmp.u+sqrt(adiabatic_term)*Pse.u));
                %Pump_combined_phase=atan2(imag(Pmp.u+sqrt(adiabatic_term)*Pse.u),real(Pmp.u+sqrt(adiabatic_term)*Pse.u));
                Pump=abs(Pmp.u).^2;
                Pump_phase=wrapToPi(angle(Pmp.u));
                Pulse=abs(Pse.u).^2;
                Pulse_phase=wrapToPi(angle(Pse.u));
            elseif Pse.component_number==2
            end
            
            if (Pse.force_pulse_time+Pse.pulse_totaltime)<=time %stage 3
                fprintf('               pulse ends\n')
                Pmp.time=time;
                Pmp.Update_u();
                if Pse.component_number==1
                    Pump=abs(Pmp.u).^2;
                    Pump_combined=Pump;
                elseif Pse.component_number==2
                    
                end
                
                Pse.clear_all();
                Pulse_add_flag=2;%enter state 3
            end
        end %pulse update end
        
        %record Movie
        if Pse.component_number==1
            
            if Pump_and_Pulse_mesh~=0 %for single component
                
                if M.re_psi_initialized==0
                    pump_fig=figure('position',[300 150 1400 900],'renderer','zbuffer');
                    M.re_psi_initialized=1;
                    M.mov_file_name_psi=Mov_file_name_psi;
                    M.currentFolder=currentFolder;
                    M.dataFolder=dataFolder;
                    
                    axis_pump=subplot(2,3,1);
                    mesh_pump=mesh(x,y,Pump);
                    axis square;
                    set(mesh_pump,'ZDataSource','Pump');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title(gca,'|Pump|^2');
                    
                    subplot(2,3,4);
                    mesh_pump_phase=mesh(x,y,Pump_phase);
                    axis square;
                    set(mesh_pump_phase,'ZDataSource','Pump_phase');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('Pump phase');
                    
                    subplot(2,3,2);
                    mesh_pulse=mesh(x,y,Pulse);
                    axis square;
                    set(mesh_pulse,'ZDataSource','Pulse');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('|Pulse|^2');
                    
                    subplot(2,3,5);
                    mesh_pulse_phase=mesh(x,y,Pulse_phase);
                    axis square;
                    set(mesh_pulse_phase,'ZDataSource','Pulse_phase');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('Pulse phase');
                    
                    axis_pump_combined=subplot(2,3,3);
                    mesh_pump_combined=mesh(x,y,Pump_combined);
                    axis square;
                    set(mesh_pump_combined,'ZDataSource','Pump_combined');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    if Pulse_add_flag==1
                        st=sprintf('Pulse added        t=%.1f',time);
                    else   st=sprintf('                   t=%.1f',time);
                    end
                    title(axis_pump_combined,{st;'|Pump + Pulse|^2'});
                    
                    subplot(2,3,6);
                    mesh_pump_combined_phase=mesh(x,y,Pump_combined_phase);
                    axis square;
                    set(mesh_pump_combined_phase,'ZDataSource','Pump_combined_phase');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('Pump + Pulse phase');
                    
                else
                    refreshdata(mesh_pump,'caller');
                    refreshdata(mesh_pump_phase,'caller');
                    refreshdata(mesh_pulse,'caller');
                    refreshdata(mesh_pulse_phase,'caller');
                    refreshdata(mesh_pump_combined,'caller');
                    refreshdata(mesh_pump_combined_phase,'caller');
                    
                    if Pulse_add_flag==1
                        st=sprintf('Pulse added        t=%.1f',time);
                    else   st=sprintf('                   t=%.1f',time);
                    end
                    title(axis_pump_combined,{st;'|Pump + Pulse|^2'});
                end %end of full Pump and Pulse mesh
                
            else %begin simplied pump_combined mesh -- Single
                
                if M.re_psi_initialized==0

                    pump_fig=figure('position',[300 300 1200 700],'renderer','zbuffer');
                    M.re_psi_initialized=1;
                    M.mov_file_name_psi=Mov_file_name_psi;
                    M.currentFolder=currentFolder;
                    M.dataFolder=dataFolder;
                    
                    axis_pump_combined=subplot(1,2,1);
                    mesh_pump_combined=mesh(x,y,Pump_combined);
                    axis square;
                    set(mesh_pump_combined,'ZDataSource','Pump_combined');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    if Pulse_add_flag==1
                        st=sprintf('Pulse added        t=%.1f',time);
                    else   st=sprintf('                   t=%.1f',time);
                    end
                    title(axis_pump_combined,{st;'|Pump + Pulse|^2'});
                    
                    subplot(1,2,2);
                    mesh_pump_combined_phase=mesh(x,y,Pump_combined_phase);
                    axis square;
                    set(mesh_pump_combined_phase,'ZDataSource','Pump_combined_phase');
                    set(gca,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('Pump + Pulse phase');
                    
                else
                    refreshdata(mesh_pump_combined,'caller');
                    refreshdata(mesh_pump_combined_phase,'caller');
                    
                    if Pulse_add_flag==1
                        st=sprintf('Pulse added        t=%.1f',time);
                    else   st=sprintf('                   t=%.1f',time);
                    end
                    title(axis_pump_combined,{st;'|Pump + Pulse|^2'});

                end
                
            end %end of one component case
            
        elseif Pse.component_number==2
            
            if isa(Pse,'LG_pulse')
                
                if M.re_psi_initialized==0 %psi windows and mesh initialization
                    
                    pump_fig=figure('position',[600 200 1000 700],'renderer','zbuffer');
                    M.re_psi_initialized=1;
                    
                    pump_L_plot=subplot(2,2,1);
                    mesh_pump_L=mesh(x,y,PumpL);
                    axis square;
                    set(mesh_pump_L,'ZDataSource','PumpL');
                    set(pump_L_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    if Pulse_add_flag==1%Pse.force_pulse_time>=time && (Pse.force_pulse_time+Pse.pulse_totaltime)<=time
                        st=sprintf('Pulse added        t=%.3f',time);
                    else   st=sprintf('                   t=%.3f',time);
                    end
                    title(pump_L_plot,{st;'|Pump_L|'});
                    
                    phase_pumpL_plot=subplot(2,2,2);
                    mesh_phase_pumpL=mesh(x,y,phase_pumpL);
                    axis square;
                    set(mesh_phase_pumpL,'ZDataSource','phase_pumpL');
                    set(phase_pumpL_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('phi Pump_L');
                    
                    PumpR_plot=subplot(2,2,3);
                    mesh_pump_R=mesh(x,y,PumpR);
                    axis square;
                    set(mesh_pump_R,'ZDataSource','PumpR');
                    set(PumpR_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    colorbar;
                    title('|Pump_R|');
                    
                    %plot phase
                    phase_pumpR_plot=subplot(2,2,4);
                    mesh_phase_pumpR=mesh(x,y,phase_pumpR);
                    axis square;
                    set(mesh_phase_pumpR,'ZDataSource','phase_pumpR');
                    set(phase_pumpR_plot,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi],'fontsize',Psi.font_size);
                    xlabel('x');
                    ylabel('y');
                    view(0,90);
                    zoom(Psi.zoom_factor);
                    caxis([-pi pi]);
                    colorbar;
                    title('phi Pump_R');
                    
                else
                    refreshdata(mesh_pump_L,'caller');
                    refreshdata(mesh_phase_pumpL,'caller');
                    refreshdata(mesh_pump_R,'caller');
                    refreshdata(mesh_phase_pumpR,'caller');
                    
                    if Pulse_add_flag==1
                        st=sprintf('Pulse added        t=%.3f',time);
                    else   st=sprintf('                   t=%.3f',time);
                    end
                    title(pump_L_plot,{st;'|Pump_L|'});
                end
                
            elseif isa(Pse,'Gaussian_pulse')
                
                
                
            end
            
        end %end of two-component initialize
        
        if M.record_movie~=0
            M.Mov_psi(pump_fig);
        else
            pause(1);
        end
        
    end
    
end % psi movie end

fprintf('\nProgram finishes.\n\n')

end

