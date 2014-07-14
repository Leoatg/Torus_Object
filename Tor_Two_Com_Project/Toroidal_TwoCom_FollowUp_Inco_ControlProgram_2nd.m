

%Control Program Object orientied

tic;
TotalTime=150;
N=600;XYmax=40;
hspace=2*XYmax/(N-1);
%dt=floor((hspace^2/pi)*1e3)*1e-3;
dt=5e-3;

%System parameters
J=0.5; GammaR_ratio=1.5; %GP parameter
m0=2;

Simulation_version='2nd';
component_number=2;

%Pump parameters
Pump0=3.5;
PL_Pum=0.4; 
l_Pum=5;
p_Pum=8;
w0_Pum=24;
Omega_Pum=0;

Full_simulation=0;% 0: FollowUP, 1: Full simulation 
GPU_calculation=0;%0:off, 1:Half GPU, 2:full GPU %present in the Mov_record class

%Pulse parameters
force_pulse_time=51;
pulse_totaltime=1000;
adiabatic_coefficient=0.5;%pulse adiabatic parameter
adabatic_shift=5;
%test adiabatic: %tt=0:0.1:60;%plot(tt,0.5*(1+tanh(0.1*(tt-11))))

Pulse0=3;
PL_Pse=0.65;
l_Pse=3;
p_Pse=6;
w0_Pse=w0_Pum;
Omega_Pse=0.25;

%Record psi_trim parameters
record_psi=0;
re_psi_begin=0;
re_psi_end=TotalTime;
re_psi_dt=0.1;

record_finial_state=1; %

%Record Movie parameters
record_movie=1;
record_mov_psi=1;
record_mov_sxyz=0;
re_mov_begin=50;
re_mov_end=TotalTime;
re_mov_dt=0.1;

font_size=12;
zoom_factor=1.2;
%----------------------------------
%add parent folder to search path
addpath(cd(cd('..')));
%----------------------------------
%Build objects
G=Grid_basic();
GPa=GP_parameter();
Pmp=LG_mode();
Pse=LG_pulse();
Psi=Psi_t();
M=Mov_record();

%initialize Grid
G.N=N;
G.XYmax=XYmax;
G.hspace=hspace;
G.dt=dt;
G.TotalTime=TotalTime;

%initialize GP paramter
GPa.J=J;

copy_Grid(Pmp,G);
copy_Grid(Pse,G);
copy_Grid(Psi,G);

copy_GP(Pmp,GPa);
copy_GP(Pse,GPa);
copy_GP(Psi,GPa);

clear G GPa
%----------------------------------
%check possible errors
if record_psi~=0 && length(re_psi_begin:re_psi_dt:re_psi_end)>20
   fprintf('\nWarning: Recorded psi might be too many (>20).\n\n');
   pause(5);
end
%----------------------------------
%initialize Pump parameter
Pmp.component_number=component_number;
Pmp.Pbar0=Pump0;
Pmp.PL=PL_Pum;
Pmp.l=l_Pum;
Pmp.p=p_Pum;
Pmp.w0=w0_Pum;
Pmp.Omega=Omega_Pum;

%initialize Pulse parameter
Pse.component_number=component_number;
Pse.Pbar0=Pulse0;
Pse.PL=PL_Pse;
Pse.l=l_Pse;
Pse.p=p_Pse;
Pse.w0=w0_Pse;
Pse.Omega=Omega_Pse;
Pse.ada_coe=adiabatic_coefficient;
Pse.ada_shift=adabatic_shift;
Pse.force_pulse_time=force_pulse_time;
Pse.pulse_totaltime=pulse_totaltime;

%initialize psi parameters
Psi.component_number=component_number;
Psi.Ini_m=m0;
Psi.record_psi=record_psi;
Psi.re_psi_begin=re_psi_begin;
Psi.re_psi_end=re_psi_end;
Psi.re_psi_dt=re_psi_dt;
Psi.font_size=font_size;
Psi.record_finial_state=record_finial_state;

%initialized record movie parameters
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
%----------------------------------------

%FollowUP Program
for Pulse0=1:0.5:2
    for Omega_Pse=0.1:0.1:0.3
        Pse.Pbar0=Pulse0;
        Pse.Omega=Omega_Pse;
        Toroidal_O_TwoCom_FollowUp_Inco_Function(Pmp,Pse,Psi,M);
    end
end

fprintf('Program finished\n ');

toc;


Clk=fix(clock);
fprintf('Finish time: %d:%d\n',Clk(4),Clk(5));






