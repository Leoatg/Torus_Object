

%% Control Program Object orientied
tic;
TotalTime=1200;
N=1024;XYmax=90;
hspace=2*XYmax/(N-1);
%dt=floor((hspace^2/pi)*1e3)*1e-3;
dt=5e-3;
damping_shift=5;%damping distance from the boundary
a_damping=0.6;

%% System parameters
GammaR=1.2; %GP parameter
fixed_int_r=1; int_r=XYmax-10;%auto set int_r to 0.95 damping
Ini_kind=2; %torodial initial state

Simulation_version='Draft_2_A';
component_number=1;

%% Pump parameters
%Pump0=2;
l_Pum=5;
p_Pum=0;
w0_Pum=28;
Omega_Pum=0;

R0=w0_Pum*sqrt(l_Pum/2)*1.1;

%% GPU calculation & Notifications
Full_simulation=1;% 0: FollowUP, 1: Full simulation
GPU_calculation=0;%0:off, 1:Half GPU, 2:full GPU %present in the Mov_record class

Send_mail=1;
Pop_up_notice=1;

%% Pulse setting
force_pulse_time=40000;
pulse_totaltime=30;
adiabatic_coefficient=0.1;
adabatic_shift=10;
%test adiabatic:%tt=0:0.1:60;%plot(tt,0.5*(1+tanh(0.1*(tt-11))))

%% LG Pulse
Pulse0=3;
l_Pse=3;
p_Pse=6;
w0_Pse=w0_Pum;
Omega_Pse=0.001;

%% Record psi_trim parameters
record_psi=0;
re_psi_begin=300;
re_psi_end=400;
re_psi_dt=5;


%% Record Movie parameters
record_movie=0;
record_mov_psi=1;
record_mov_sxyz=0;
re_mov_begin=300;
re_mov_end=TotalTime;
re_mov_dt=0.1;

font_size=12;
zoom_factor=1;

%% add parent folder to search path
addpath(cd(cd('..')));
addpath('D:\Torus_Pulse_Object');

%% Main parameter loop
for Pump0=2:-0.1:2
%for GammaR=1.2:0.2:1.2
for m0=40:40
    
%% Build objects
G=Grid_basic();
GPa=GP_parameter();
Pmp=LG_mode();
Pse=LG_pulse();
Psi=Psi_t();
M=Mov_record();

%% initialize Grid
G.N=N;
G.XYmax=XYmax;
G.hspace=hspace;
G.dt=dt;
G.TotalTime=TotalTime;
G.damping_shift=damping_shift;%damping distance from the boundary
G.a_damping=a_damping;

%% initialize GP paramter
GPa.GammaR=GammaR;
GPa.GammaR_ratio=GammaR/GPa.GammaC;
GPa.Update_GP_parameter();

copy_Grid(Pmp,G);
copy_Grid(Pse,G);
copy_Grid(Psi,G);

copy_GP(Pmp,GPa);
copy_GP(Pse,GPa);
copy_GP(Psi,GPa);

clear G GPa

%% check possible errors
if record_psi~=0 && length(re_psi_begin:re_psi_dt:re_psi_end)>20
    fprintf('\nWarning: Recorded psi might be too many (>20).\n\n');
    pause(0.5);
end

%% initialize Pump parameter
Pmp.component_number=component_number;
Pmp.Pbar0=Pump0;
Pmp.l=l_Pum;
Pmp.p=p_Pum;
Pmp.w0=w0_Pum;
Pmp.Omega=Omega_Pum;

%% initialize Pulse parameter
Pse.component_number=component_number;
Pse.Pbar0=Pulse0;
Pse.l=l_Pse;
Pse.p=p_Pse;
Pse.w0=w0_Pse;
Pse.Omega=Omega_Pse;
Pse.ada_coe=adiabatic_coefficient;
Pse.ada_shift=adabatic_shift;
Pse.force_pulse_time=force_pulse_time;
Pse.pulse_totaltime=pulse_totaltime;

%% initialize psi parameters
Psi.component_number=component_number;
Psi.Ini_m=m0;
Psi.fixed_int_r=fixed_int_r;
Psi.int_r=int_r;
Psi.Ini_R=R0;
Psi.record_psi=record_psi;
Psi.re_psi_begin=re_psi_begin;
Psi.re_psi_end=re_psi_end;
Psi.re_psi_dt=re_psi_dt;
Psi.font_size=font_size;

%% initialized record movie parameters
M.record_movie=record_movie;
M.re_mov_begin=re_mov_begin;
M.re_mov_end=re_mov_end;
M.re_mov_dt=re_mov_dt;
M.record_mov_psi=record_mov_psi;
M.record_mov_sxyz=record_mov_sxyz;
M.zoom_factor=zoom_factor;
Psi.zoom_factor=zoom_factor;
M.GPU_calculation=GPU_calculation;
M.Simulation_version=Simulation_version;

%% Full Simulation

Toroidal_O_Single_Full_Simulation_Inco_Function(Pmp,Pse,Psi,M);

end
end
%% Notifications
%----------------------------------------
Clk=fix(clock);
st_finish_time=sprintf('Finish time: %d.%d %d:%d \n',Clk(2),Clk(3),Clk(4),Clk(5));
%----------------------------------------
if Send_mail~=0 || Pop_up_notice~=0
    
    st_mail=['Program: ',Simulation_version,sprintf('_T=%d',floor(TotalTime))];
    
    if record_movie~=0
        st_mail=[st_mail,'_Mov'];
    end
    if GPU_calculation~=0
        st_mail=[st_mail,'_GPU'];
    end
    
    st_mail={[st_mail,' finished.'],st_finish_time};
    
    if Send_mail~=0
        %send mail code
        setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');
        setpref('Internet','E_mail','matlab_notification@yahoo.com.au');
        props = java.lang.System.getProperties;
        props.setProperty('mail.smtp.auth','true');
        setpref('Internet','SMTP_Username','matlab_notification@yahoo.com.au');
        setpref('Internet','SMTP_Password','Aaa12345');
        props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
        props.setProperty('mail.smtp.socketFactory.port', '465');
        sendmail('chinaliguang@msn.com','Automatic Generated Matlab Email',st_mail);
    end
end
%----------------------------------------
%pop up message
if Pop_up_notice~=0
    h_msbox = msgbox(st_mail,'Done');
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jf=get(get(h_msbox,'JavaFrame'),'FigurePanelContainer');%set always on top
    jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
end
%----------------------------------------

fprintf('Program finished\n ');

toc;

fprintf(st_finish_time);








