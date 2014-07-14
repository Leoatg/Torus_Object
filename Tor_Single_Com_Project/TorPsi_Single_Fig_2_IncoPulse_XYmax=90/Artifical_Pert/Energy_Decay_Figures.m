

%energy decay figures

addpath('d:\Torus_Pulse_Object');

if ~exist('psi_trim','var')
load('Tor_Pert_GmaRr=1.4_m0=76_Pmp0=2.5_l=5_p=0_w0=28_PerS=0.05_PerF=-56_PerT=60_XYmax=90_N=2048_T=301.mat');
end

%build x y
XYmax=Ori_Pump_Parameter.XYmax;
N=Ori_Pump_Parameter.N;
hspace=2*XYmax/(N-1);
[x,y]=meshgrid(-XYmax:hspace:XYmax,-XYmax:hspace:XYmax);

Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;
Windows_shift=100;
zoom_factor=1;

print_fig=1;

%plot range
x_b=70;
N_b=floor((XYmax-x_b)/hspace);
N_e=floor((XYmax+x_b)/hspace);

x_t=x(N_b:N_e,N_b:N_e);
y_t=y(N_b:N_e,N_b:N_e);

%t=300, 340, 400
%   1    9    21
psi1=psi_trim(:,:,1);
psi2=psi_trim(:,:,13);
psi3=psi_trim(:,:,19);

fig_1_den=figure('position',[362   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');

mesh(x_t,y_t,abs(psi1(N_b:N_e,N_b:N_e,1)).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(c)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


pos_1=get(fig_1_den,'position');

fig_1_phase=figure('position',[pos_1(1)+Windows_shift   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,wrapToPi(angle(psi1(N_b:N_e,N_b:N_e,1))));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'zlim',[-pi,pi]);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(d)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


fig_2_den=figure('position',[pos_1(1)+Windows_shift*2   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,abs(psi2(N_b:N_e,N_b:N_e,1)).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(e)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


fig_2_phase=figure('position',[pos_1(1)+Windows_shift*3   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,wrapToPi(angle(psi2(N_b:N_e,N_b:N_e,1))));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'zlim',[-pi,pi]);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(f)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


fig_3_den=figure('position',[pos_1(1)+Windows_shift*4   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,abs(psi3(N_b:N_e,N_b:N_e,1)).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(g)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


fig_3_phase=figure('position',[pos_1(1)+Windows_shift*5   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,wrapToPi(angle(psi3(N_b:N_e,N_b:N_e,1))));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'zlim',[-pi,pi]);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+12;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(h)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


fig_E=figure('position',[pos_1(1)+Windows_shift*6   333   windows_width+200   windows_hight+100],'renderer','painter','paperpositionmode','auto');
[AX,H1,H2]=plotyy(Time,Energy,Time,real(Lz));
set(AX(1),'fontsize',Font_size,'fontname',Font_name);
set(AX(1),'xlim',[0,max(Time)]);
set(AX(1),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(1),'ylim',[1 7],'ytick',1:7);
set(H1,'linewidth',3);
xlabel('t');
ylabel(AX(1),'Energy');
box(AX(1),'off');

set(AX(2),'fontsize',Font_size,'fontname',Font_name);
set(AX(2),'xlim',[0,max(Time)]);
set(AX(2),'xtick',[]);
set(AX(2),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'ylim',[10 80]);
set(AX(2),'ytick',10:20:80);
set(H2,'linewidth',3,'linestyle','--');
ylabel(AX(2),'L_z');
box(AX(2),'off');

linkaxes(AX,'x');
set(AX(2), 'XTickLabel','','XAxisLocation','Top') 

ann_position=[ 0.7531    0.8267    0.0800    0.0800];
h_a=annotation('textbox', ann_position,'String','(b)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+4,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',ann_position);

if print_fig~=0
    
   resolution='-r180';
   print(fig_1_den,'-zbuffer','-depsc2',resolution,'Energy_Decay_Den_1');
   print(fig_2_den,'-zbuffer','-depsc2',resolution,'Energy_Decay_Den_2');
   print(fig_3_den,'-zbuffer','-depsc2',resolution,'Energy_Decay_Den_3');
   print(fig_1_phase,'-zbuffer','-depsc2',resolution,'Energy_Decay_phase_1');
   print(fig_2_phase,'-zbuffer','-depsc2',resolution,'Energy_Decay_phase_2');
   print(fig_3_phase,'-zbuffer','-depsc2',resolution,'Energy_Decay_phase_3');
   
   print(fig_E,'-painter','-depsc','-r250','Energy_Decay_Energy.eps');
   
   close all;    
end














