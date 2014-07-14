

%Control Program Object orientied
%Full simulation 
tic;
TotalTime=130;
N=1500;XYmax=35;
hspace=2*XYmax/(N-1);
%dt=floor((hspace^2/pi)*1e4)*1e-4;
dt=2e-3;
damping_shift=5;%damping distance from the boundary
a_damping=0.5; 

%System parameters
J=0.5; GammaR_ratio=1.5; %GP parameter
Init_ratio=0.1; %Psi_t parameter
m0=0;%Initial angular momentum
Ini_kind=1;
fixed_int_r=1; int_r=18;
int_tolerant=1e-3;

Simulation_version='1st_Elo';
component_number=2;

%Pump parameters
Pump0=3;
PL_Pmp=0.4;
Sigma_x_Pmp=16;
Sigma_y_Pmp=5;
x0_Pmp=0;
y0_Pmp=0;

GPU_calculation=0;%0:off, 1:Half GPU, 2:full GPU %present in the Mov_record class

%Pulse parameters
force_pulse_time=800;
pulse_totaltime=100;
adiabatic_coefficient=0.5;%pulse adiabatic parameter
adabatic_shift=5;%test adiabatic: %tt=0:0.1:60;%plot(tt,0.5*(1+tanh(0.1*(tt-11))))

Pulse0=3;
PL_Pse=0.65;
Sigma_x_Pse=5;
Sigma_y_Pse=5;
x0_Pse=0;
y0_Pse=0;

%Record psi_trim parameters
record_psi=0;
re_psi_begin=0;
re_psi_end=TotalTime;
re_psi_dt=0.1;

record_finial_state=1; %
record_sxyz_int=0;

%Record Movie parameters
record_movie=0;
record_mov_psi=1;
record_mov_sxyz=1;
re_mov_begin=0;
re_mov_end=TotalTime;
re_mov_dt=0.1;

font_size=12;
zoom_factor=1.2;
%----------------------------------
%add parent folder to search path
addpath(cd(cd('..')));
addpath('D:\Torus_Pulse_Object');
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
G.damping_shift=damping_shift;%damping distance from the boundary
G.a_damping=a_damping; 

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
    %return;
elseif (Sigma_x_Pmp==Sigma_y_Pmp) && Pump0<4
    fprintf('\nError: Pumping value is too weak for the round pump\n\n');
    %return;
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
Psi.Ini_kind=Ini_kind;
Psi.fixed_int_r=fixed_int_r;
Psi.Ini_m=m0;
Psi.fixed_int_r=fixed_int_r;
Psi.int_r=int_r;
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

Font_size=16;
zoom_factor=1;

PhaseSlop=x+2*y;
Px0=0;
Py0=1;
Phase_V=atan2(y-Py0,x-Px0);

figure('renderer','zbuffer','position',[50 100 1800 1000],'paperpositionmode','auto');

subplot(2,4,1);
mesh(x,y,Phase_V);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
title('original vortex');

subplot(2,4,5);
mesh(x,y,wrapToPi(Phase_V));
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
view(0,90);
zoom(zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('wrapped  vortex');

subplot(2,4,2);
mesh(x,y,PhaseSlop);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
title('original phase slop');

subplot(2,4,3);
mesh(x,y,wrapToPi(PhaseSlop));
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
view(0,90);
zoom(zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('wrapped original phase slop');

subplot(2,4,4);
mesh(x,y,unwrap(unwrap(wrapToPi(PhaseSlop),[],1),[],2));
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
title('unwrap the wrapped phase slop by Matlab');


subplot(2,4,6);
mesh(x,y,PhaseSlop+Phase_V);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
title('unwrapped phase slop with vortex');

subplot(2,4,7);
mesh(x,y,wrapToPi(PhaseSlop+Phase_V));
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
view(0,90);
zoom(zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('wrapped phase slop with vortex');

subplot(2,4,8);
mesh(x,y,unwrap(unwrap(wrapToPi(PhaseSlop+Phase_V),[],1),[],2));
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
title('unwrap the wrapped phase slop with vortex by Matlab');

print(gcf,'-zbuffer','-dpng','-r250','phase_slop_and_vortex.png');


fprintf('Program finished\n ');

toc;







