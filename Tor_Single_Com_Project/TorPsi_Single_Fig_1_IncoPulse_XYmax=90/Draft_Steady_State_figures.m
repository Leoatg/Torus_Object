
%Draft 1 steady state figures

addpath('d:\Torus_Pulse_Object');

%load('FinalState_TorPsi_IncoPulse_GmaRr=1.5_J=0.5_m0=0.0_Pmp0=2.50_l=5_p=0_w0=28_XYmax=90_N=2048_T=80.mat');
load('FinalState_TorPsi_IncoPulse_GmaRr=1.5_m0=0.0_Pmp0=2.50_l=5_p=0_w0=28_XYmax=90_N=2048_T=81.mat');

%build x y
XYmax=Ori_Pump_Parameter.XYmax;
N=Ori_Pump_Parameter.N;
hspace=2*XYmax/(N-1);
[x,y]=meshgrid(-XYmax:hspace:XYmax,-XYmax:hspace:XYmax);
x_axis=-XYmax:hspace:XYmax;

%build pump
Pmp=LG_mode();
Pmp.XYmax=XYmax;
Pmp.N=N;
Pmp.hspace=hspace;
Pmp.Pbar0=Ori_Pump_Parameter.Pbar0;
Pmp.l=Ori_Pump_Parameter.l;
Pmp.p=Ori_Pump_Parameter.p;
Pmp.w0=Ori_Pump_Parameter.w0;
Pmp.Build_u();
Pump=abs(Pmp.u).^2;


Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;
Windows_shift=100;
zoom_factor=1.1;

%plot range
x_b=70;
N_b=floor((XYmax-x_b)/hspace);
N_e=floor((XYmax+x_b)/hspace);

x_t=x(N_b:N_e,N_b:N_e);
y_t=y(N_b:N_e,N_b:N_e);

print_fig=1;


fig_1b=figure('position',[362   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');

%mesh(x,y,abs(Finial_State(:,:,1)).^2);
mesh(x_t,y_t,abs(Finial_State(N_b:N_e,N_b:N_e,1)).^2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
%set(gca,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
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

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(b)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);

pos_1=get(fig_1b,'position');

fig_1c=figure('position',[pos_1(1)+Windows_shift   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
%mesh(x,y,wrapToPi(angle(Finial_State(:,:,1))));
mesh(x_t,y_t,wrapToPi(angle(Finial_State(N_b:N_e,N_b:N_e,1))));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
%set(gca,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax],'zlim',[-pi,pi]);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b],'zlim',[-pi,pi]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+4;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(c)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);

%radial profile & Pump
fig_radial_pro=figure('position',[pos_1(1)+Windows_shift*3   333   windows_width   windows_hight],'renderer','painter','paperpositionmode','auto');
%[AX,H1,H2]=plotyy(x_axis,Pump(N/2+1,:),x_axis,abs(Finial_State(N/2+1,:,1)).^2);
[AX,H1,H2]=plotyy(x_axis(N/2+1:end),Pump(N/2+1,N/2+1:end),x_axis(N/2+1:end),abs(Finial_State(N/2+1,N/2+1:end,1)).^2);
axis square;
set(AX(1),'fontsize',Font_size,'fontname',Font_name);
set(AX(2),'fontsize',Font_size,'fontname',Font_name,'xtick',[],'PlotBoxAspectRatio',[1,1,1]);
set(H1,'linewidth',2.5);
set(H2,'linewidth',2.5,'linestyle','--')
xlabel(AX(1),'x');
ylabel(AX(1),'Pump');
ylabel(AX(2),'$\Phi^2$','interpreter','latex');



load('FinalState_TorPsi_IncoPulse_GmaRr=1.5_J=0.5_m0=2.0_Pmp0=2.50_l=5_p=0_w0=28_XYmax=90_N=2048_T=80.mat','Finial_State');

fig_1d=figure('position',[pos_1(1)+Windows_shift*2   333   windows_width   windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(x_t,y_t,wrapToPi(angle(Finial_State(N_b:N_e,N_b:N_e,1))));
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
%set(gca,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax],'zlim',[-pi,pi]);
set(gca,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b],'zlim',[-pi,pi]);
view(0,90);
zoom(zoom_factor);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
xlabel('x');
h_y=ylabel('y');
pos_y=get(h_y,'position');
pos_y(1)=pos_y(1)+4;
set(h_y,'position',pos_y);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(d)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);

if print_fig~=0
    
    currentFolder = pwd;
    dataFolder_name='Draft_Fig_Standard_Steady_State';
    
    if exist(dataFolder_name,'dir')==0
        mkdir(dataFolder_name);
    end
    
    
    resolution='-r180';
    print(fig_1b,'-zbuffer','-depsc2',resolution,[currentFolder,'\',dataFolder_name,'\','Steady_state_mode_m0.eps']);
    print(fig_1c,'-zbuffer','-depsc2',resolution,[currentFolder,'\',dataFolder_name,'\','Steady_state_phase_m0.eps']);
    print(fig_1d,'-zbuffer','-depsc2',resolution,[currentFolder,'\',dataFolder_name,'\','Steady_state_phase_m2.eps']);
    
    print(fig_radial_pro,'-painter','-depsc2','-r300',[currentFolder,'\',dataFolder_name,'\','Steady_state_Density_Profile.eps']);
    
    
    close all;
end

















