

%Control Program Object orientied

tic;
TotalTime=5;
N=1024;XYmax=90;
hspace=2*XYmax/(N-1);
dt=floor((hspace^2/pi)*1e3)*1e-3;
%dt=1;
%System parameters
J=0.5; GammaR_ratio=1.5; %GP parameter
Init_ratio=0.1; %Psi_t parameter

Simulation_version='test';
component_number=1;

%Pump parameters
Pump0=5;
PL_Pum=0.4; 
l_Pum=5;
p_Pum=0;
w0_Pum=8;
Omega_Pum=0;
z_Pum=40;%zR=5

Pulse_kind=0; %0: LG pulse, 1: Gaussian pulse

%Pulse parameters
force_pulse_time=0;
pulse_totaltime=1000;
adiabatic_coefficient=0.5;%pulse adiabatic parameter
adabatic_shift=5;

%LG Pulse parameter
Pulse0_LG=1;
PL_Pse=0.65;
l_Pse=-5;
p_Pse=0;
w0_Pse=16;
Omega_LG=2e-2;
z_Pse=z_Pum;

%Gaussian Pulse parameter
Pulse0_Gau=2;
x0=20;
y0=0;
Omega_Gau=0.1;
Sigma_x=5;
Sigma_y=5;

%Record Movie parameters
record_movie=0;
record_mov_psi=1;
record_mov_sxyz=0;
re_mov_begin=0;
re_mov_end=TotalTime;
re_mov_dt=10;

font_size=12;
zoom_factor=1;
%----------------------------------
%Build objects
G=Grid_basic();
GPa=GP_parameter();
Pmp=LG_mode();
if Pulse_kind==0
    Pse=LG_pulse();
elseif Pulse_kind==1
    Pse=Gaussian_pulse();
end
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

%initialize Pump parameter
Pmp.component_number=component_number;
Pmp.Pbar0=Pump0;
Pmp.PL=PL_Pum;
Pmp.l=l_Pum;
Pmp.p=p_Pum;
Pmp.w0=w0_Pum;
Pmp.Omega=Omega_Pum;
Pmp.z=z_Pum;

%initialize Pulse parameter
Pse.component_number=component_number;
Pse.force_pulse_time=force_pulse_time;
Pse.pulse_totaltime=pulse_totaltime;
Pse.ada_coe=adiabatic_coefficient;
Pse.ada_shift=adabatic_shift;
if Pulse_kind==0
    Pse.Pbar0=Pulse0_LG;
    Pse.l=l_Pse;
    Pse.p=p_Pse;
    Pse.w0=w0_Pse;
    Pse.Omega=Omega_LG;
    Pse.z=z_Pse;
    if component_number==2
        Pse.PL=PL_Pse;
    end
elseif Pulse_kind==1
    Pse.Pbar0=Pulse0_Gau;
    Pse.Omega=Omega_Gau;
    Pse.x0=x0;
    Pse.y0=y0;
    Pse.Sigma_x=Sigma_x;
    Pse.Sigma_y=Sigma_y;
    if component_number==2
        Pse.PL=PL_Pse;
    end
end

%initialized record movie parameters
M.record_movie=record_movie;
M.re_mov_begin=re_mov_begin;
M.re_mov_end=re_mov_end;
M.re_mov_dt=re_mov_dt;
M.record_mov_psi=record_mov_psi;
M.record_mov_sxyz=record_mov_sxyz;
M.zoom_factor=zoom_factor;
Psi.zoom_factor=zoom_factor;
M.Simulation_version=Simulation_version;
%----------------------------------------

Toroidal_TestPump_Function(Pmp,Pse,Psi,M);

fprintf('Program finished\n ');

toc;







