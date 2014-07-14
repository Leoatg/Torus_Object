
% m_es phase diagram

ua=7.7e-3;
R1=8.4e-3;

P1=1:0.001:5;
GammaR=0.1:0.001:4.002;

[x,y]=meshgrid(P1,GammaR);

Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;

print_fig=1;


fig_mes_diag=figure('position',[300 200 windows_width windows_hight],'renderer','zbuffer','paperpositionmode','auto');

mesh(x,y,sqrt(ua*y.*(x-1)/R1));
axis square;
set(gca,'xlim',[min(P1),max(P1)],'ylim',[min(GammaR),max(GammaR)],'fontsize',Font_size,'Fontname',Font_name);
view(0,90);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
h_x=xlabel('$\bar{P}$','Interpreter','LaTex');
set(h_x,'position',[3.3   -0.2    2.0001]);
ylabel('$\gamma_R$','Interpreter','Latex');

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(a)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.750    0.08    0.08]);



if print_fig~=0
    
    resolution='-r180';
    
    print(fig_mes_diag,'-zbuffer','-depsc2',resolution,'m_es_Phase_Diagram.eps');
    print(fig_mes_diag,'-zbuffer','-dpng',resolution,'m_es_Phase_Diagram.png');
    
    close(fig_mes_diag);
end





















