
%Plot m_ds Phase Diagram

load('m_ds_phase_Diagram_Fine_4_R0=44.mat');

Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;

print_fig=1;


[x,y]=meshgrid(P1_series,GaR_series);

%interp
[X,Y]=meshgrid(P1_min:0.001:P1_max,GaR_min:0.001:GaR_max);
[m_ds_inted]=interp2(x,y,m_ds,X,Y);

fig_mds_diag2=figure('position',[300 200 windows_width windows_hight],'renderer','zbuffer','paperpositionmode','auto');
mesh(X,Y,m_ds_inted/R0);
axis square;
set(gca,'xlim',[P1_min,P1_max],'ylim',[GaR_min,GaR_max]);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
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
    
    % print(fig_mes_diag,'-zbuffer','-depsc',resolution,'m_ds_Phase_Diagram.eps');
    %print(fig_mds_diag,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram.png');
    %print(fig_mds_diag2,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram2.png');
    print(fig_mds_diag2,'-zbuffer','-depsc2',resolution,sprintf('m_ds_Phase_Diagram_R0=%.1f.eps',R0));
    
    close(fig_mes_diag2);
end


fprintf('Program finished.\n');

toc;