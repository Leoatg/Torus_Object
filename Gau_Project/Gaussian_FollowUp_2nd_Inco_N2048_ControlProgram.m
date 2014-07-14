

%Control Program Object orientied
%Full simulation 
tic;
TotalTime=180;
N=2048;XYmax=60;
hspace=2*XYmax/(N-1);
dt=floor((hspace^2/pi)*1e3)*1e-3;

%System parameters
J=0.5; GammaR_ratio=1.5; %GP parameter
Init_ratio=0.1; %Psi_t parameter
m0=0;%Initial angular momentum
int_tolerant=1e-3;

Simulation_version='2nd_N2048';

%Pump parameters
Pump0=3;
PL_Pmp=0.4;
Sigma_x_Pmp=17;
Sigma_y_Pmp=5;
x0_Pmp=0;
y0_Pmp=0;

GPU_calculation=0;%0:off, 1:Half GPU, 2:full GPU %present in the Mov_record class

%Pulse parameters
force_pulse_time=130;
pulse_totaltime=100;
adiabatic_coefficient=0.5;%pulse adiabatic parameter
adabatic_shift=5;%test adiabatic: %tt=0:0.1:60;%plot(tt,0.5*(1+tanh(0.1*(tt-11))))

Pulse0=3;
PL_Pse=0.65;
Sigma_x_Pse=4;
Sigma_y_Pse=4;
x0_Pse=0;
y0_Pse=0;

%Record psi_trim parameters
record_psi=1;
re_psi_begin=170;
re_psi_end=180;
re_psi_dt=0.5;

record_finial_state=0; %
record_sxyz_int=0;

%Record Movie parameters
record_movie=0;
record_mov_psi=1;
record_mov_sxyz=1;
re_mov_begin=0;
re_mov_end=TotalTime;
re_mov_dt=0.5;

font_size=12;
zoom_factor=1.5;
%----------------------------------
%add parent folder to search path
addpath(cd(cd('..')));
%----------------------------------
%Build objects
G=Grid_basic();
GPa=GP_parameter();
Pmp=Gaussian_mode();
Pse=Gaussian_pulse();
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
if (Sigma_x_Pmp~=Sigma_y_Pmp) && Pump0>3
    fprintf('\nError: Elongated pumping value beyond the synchronized threshold\n\n');
    return;
elseif (Sigma_x_Pmp==Sigma_y_Pmp) && Pump0<4
    fprintf('\nError: Pumping value is too weak for the round pump\n\n');
    return;
end
if record_psi~=0 && length(re_psi_begin:re_psi_dt:re_psi_end)>20
   fprintf('\nWarning: Recorded psi might be too many (>20).\n\n');
   pause(5);
end
%set integration area
if Sigma_x_Pmp==Sigma_y_Pmp
   Psi.fixed_int_r=0;
   Psi.int_tolerant=int_tolerant;
end
%----------------------------------
%initialize Pump parameter
Pmp.Pbar0=Pump0;
Pmp.PL=PL_Pmp;
Pmp.Sigma_x=Sigma_x_Pmp;
Pmp.Sigma_y=Sigma_y_Pmp;
Pmp.x0=x0_Pmp;
Pmp.y0=y0_Pmp;

%initialize Pulse parameter
Pse.Pbar0=Pulse0;
Pse.PL=PL_Pse;
Pse.Sigma_x=Sigma_x_Pse;
Pse.Sigma_y=Sigma_y_Pse;
Pse.x0=x0_Pse;
Pse.y0=y0_Pse;
Pse.ada_coe=adiabatic_coefficient;
Pse.ada_shift=adabatic_shift;
Pse.force_pulse_time=force_pulse_time;
Pse.pulse_totaltime=pulse_totaltime;

%initialize psi parameters
Psi.IR=Init_ratio;
Psi.Ini_m=m0;
Psi.record_psi=record_psi;
Psi.re_psi_begin=re_psi_begin;
Psi.re_psi_end=re_psi_end;
Psi.re_psi_dt=re_psi_dt;
Psi.record_sxyz_int=record_sxyz_int;
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

%Full simulation

Gaussian_O_FollowUp_Inco_Function(Pmp,Pse,Psi,M);

fprintf('Program finished\n ');

toc;







