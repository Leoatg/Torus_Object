

%energy decay figures

addpath('d:\Torus_Pulse_Object');

if ~exist('psi_trim','var')
load('TorPsi_IncoPulse_GmaRr=1.5_m0=44.0_Pmp0=2.50_l=5_p=0_w0=16_XYmax=90_N=2048_T=400.mat');
end

%build x y
XYmax=Ori_Pump_Parameter.XYmax;
N=Ori_Pump_Parameter.N;
hspace=2*XYmax/(N-1);
[x,y]=meshgrid(-XYmax:hspace:XYmax,-XYmax:hspace:XYmax);

Font_size=22;
Font_name='Times New Roman';
Windows_shift=25;
zoom_factor=1.1;

print_fig=0;

%t=300, 340, 400
%   1    9    21
psi1=psi_trim(:,:,1);
psi2=psi_trim(:,:,9);
psi3=psi_trim(:,:,21);

fig_1_den=figure('position',[362   333   805   649],'renderer','zbuffer','paperpositionmode','auto');

mesh(x,y,abs(psi1).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

pos_1=get(fig_1_den,'position');

fig_1_phase=figure('position',[pos_1(1)+Windows_shift   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
mesh(x,y,wrapToPi(angle(psi1)));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax],'zlim',[-pi,pi]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

fig_2_den=figure('position',[pos_1(1)+Windows_shift*2   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
mesh(x,y,abs(psi2).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

fig_2_phase=figure('position',[pos_1(1)+Windows_shift*3   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
mesh(x,y,wrapToPi(angle(psi2)));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax],'zlim',[-pi,pi]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

fig_3_den=figure('position',[pos_1(1)+Windows_shift*4   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
mesh(x,y,abs(psi3).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

fig_3_phase=figure('position',[pos_1(1)+Windows_shift*5   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
mesh(x,y,wrapToPi(angle(psi3)));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax],'zlim',[-pi,pi]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);


fig_E=figure('position',[pos_1(1)+Windows_shift*6   333   805   649],'renderer','painter','paperpositionmode','auto');
[AX,H1,H2]=plotyy(Time,Energy,Time,real(Lz));
set(AX(1),'fontsize',Font_size,'fontname',Font_name);
set(AX(1),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(1),'ylim',[1 7]);
set(H1,'linewidth',3);
xlabel('t');
ylabel(AX(1),'Energy');

set(AX(2),'fontsize',Font_size,'fontname',Font_name);
set(AX(2),'xtick',[]);
set(AX(2),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'ylim',[-10 55]);
set(H2,'linewidth',3,'linestyle','--');
ylabel(AX(2),'L_z');

if print_fig~=0
    
   resolution='-r100';
   print(fig_1_den,'-zbuffer','-depsc',resolution,'Energy_Decay_Den_1');
   print(fig_2_den,'-zbuffer','-depsc',resolution,'Energy_Decay_Den_2');
   print(fig_3_den,'-zbuffer','-depsc',resolution,'Energy_Decay_Den_3');
   print(fig_1_phase,'-zbuffer','-depsc',resolution,'Energy_Decay_phase_1');
   print(fig_2_phase,'-zbuffer','-depsc',resolution,'Energy_Decay_phase_2');
   print(fig_3_phase,'-zbuffer','-depsc',resolution,'Energy_Decay_phase_3');
   
   print(fig_E,'-painter','-depsc','-r250','Energy_Decay_Energy.eps');
   
   close all;    
end














