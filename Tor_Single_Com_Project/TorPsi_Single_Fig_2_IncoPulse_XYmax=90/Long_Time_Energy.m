
% long time Energy

addpath('d:\Torus_Pulse_Object');

load('FinalState_TorPsi_IncoPulse_GmaRr=1.5_m0=44.0_Pmp0=2.50_l=5_p=0_w0=16_XYmax=90_N=1024_T=500.mat');

%build x y
XYmax=Ori_Pump_Parameter.XYmax;
N=Ori_Pump_Parameter.N;
hspace=2*XYmax/(N-1);
[x,y]=meshgrid(-XYmax:hspace:XYmax,-XYmax:hspace:XYmax);

Font_size=22;
Font_name='Times New Roman';
Windows_shift=25;
zoom_factor=1.1;

print_fig=1;

fig_E=figure('position',[362   333   805   649],'renderer','painter','paperpositionmode','auto');
[AX,H1,H2]=plotyy(Time,Energy,Time,real(Lz));
set(AX(1),'fontsize',Font_size,'fontname',Font_name);
set(AX(1),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(H1,'linewidth',3);
xlabel('t');
ylabel(AX(1),'Energy');

set(AX(2),'fontsize',Font_size,'fontname',Font_name);
set(AX(2),'xtick',[]);
set(AX(2),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(H2,'linewidth',3,'linestyle','--');
ylabel(AX(2),'L_z');


if print_fig~=0

print(fig_E,'-painter','-depsc','-r250','Energy_Decay_Energy.eps');
   
close(fig_E);




end


