

%energy decay figures

addpath('d:\Torus_Pulse_Object');

if ~exist('psi_trim','var')
load('FinalState_TorPsi_IncoPulse_GmaRr=0.5_m0=0.0_Pmp0=2.50_l=5_p=0_w0=28_XYmax=90_N=1024_T=240.mat');
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

%t=240
%   10   
psi=Finial_State(:,:,1);
n=Finial_State(:,:,2);

phi=wrapToPi(angle(psi));

%plot finial density
fig_den=figure('position',[362   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');

mesh(x_t,y_t,abs(psi(N_b:N_e,N_b:N_e,1)).^2);
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
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(c)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);


pos_1=get(fig_den,'position');

%plot finial phase
fig_phase=figure('position',[pos_1(1)+Windows_shift   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,phi(N_b:N_e,N_b:N_e,1));
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
pos_y(1)=pos_y(1)+8;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(d)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);

%plot reservoir
% fig_n=figure('position',[pos_1(1)+Windows_shift*2   333   805   649],'renderer','zbuffer','paperpositionmode','auto');
% mesh(x_t,y_t,n(N_b:N_e,N_b:N_e,1));
% axis square;
% set(gca,'fontsize',Font_size,'fontname',Font_name);
% set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
% view(0,90);
% zoom(zoom_factor);
% colormap(jet(4096));
% shading interp;
% colorbar('fontsize',Font_size,'fontname',Font_name);
% xlabel('x');
% h_y=ylabel('y');
% pos_y=get(h_y,'position');
% pos_y(1)=pos_y(1)+8;
% set(h_y,'position',pos_y);

%plot energy
% fig_E=figure('position',[pos_1(1)+Windows_shift*3   333   805   649],'renderer','painter','paperpositionmode','auto');
% plot(Time,Energy,'linewidth',3);
% axis tight;
% set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[min(Time),max(Time)]);
% xlabel('t');
% ylabel('Energy');

%plot Lz
fig_Lz=figure('position',[pos_1(1)+Windows_shift*3   333   windows_width-50   windows_hight-50],'renderer','painter','paperpositionmode','auto');
plot(Time,real(Lz),'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name,'xlim',[0,max(Time)],'ylim',[-1.5,2.5]);
xlabel('t');
ylabel('Lz');

ann_position=[0.1920    0.8084    0.0800    0.0800];
h_a=annotation('textbox', ann_position,'String','(b)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',ann_position);


% fig_E=figure('position',[pos_1(1)+Windows_shift*6   333   805   649],'renderer','painter','paperpositionmode','auto');
% [AX,H1,H2]=plotyy(Time,Energy,Time,real(Lz));
% set(AX(1),'fontsize',Font_size,'fontname',Font_name);
% set(AX(1),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
% set(H1,'linewidth',3);
% xlabel('t');
% ylabel(AX(1),'Energy');
% 
% set(AX(2),'fontsize',Font_size,'fontname',Font_name);
% set(AX(2),'xtick',[]);
% set(AX(2),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
% set(H2,'linewidth',3,'linestyle','--');
% ylabel(AX(2),'L_z');



if print_fig~=0
    
   resolution='-r180';
   print(fig_den,'-zbuffer','-depsc2',resolution,'Dynamic_Decay_Den');
   %print(fig_n,'-zbuffer','-depsc2',resolution,'Dynamic_Decay_n');
   print(fig_phase,'-zbuffer','-depsc2',resolution,'Dynamic_Decay_phase');
     
   %print(fig_E,'-painter','-depsc2','-r250','Dynamic_Decay_Energy.eps');
   print(fig_Lz,'-painter','-depsc2','-r300','Dynamic_Decay_Lz.eps');
   
   close all;    
end














