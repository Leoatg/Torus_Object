
% m_es phase diagram

ua=7.7e-3;
R1=8.4e-3;

P1=2.107;
GammaR=1.4;

Phi0=sqrt(GammaR*(P1-1)/R1);

k=-80:0.01:20;

Font_size=20;
Font_name='Times New Roman';

print_fig=0;

R0=24*sqrt(5/2);

m=64;

omega=(2*k*m+sqrt(k.^4+k.^2*4*ua*R0^2*Phi0^2))/(2*R0^2);

figure('renderer','painter','paperpositionmode','auto');
plot(k,omega,'linewidth',2);
axis square;
xlabel('k');
ylabel('\omega');
hold all
plot(k,zeros(1,length(k)));
hold off

[a,b]=min(omega);

min_k=k(b)


% fig_mes_diag=figure('position',[300 200 800 600],'renderer','zbuffer','paperpositionmode','auto');
% 
% mesh(x,y,sqrt(ua*y.*(x-1)/R1));
% axis square;
% set(gca,'xlim',[min(P1),max(P1)],'ylim',[min(GammaR),max(GammaR)],'fontsize',Font_size,'Fontname',Font_name);
% view(0,90);
% colormap(jet(4096));
% shading interp;
% colorbar('fontsize',Font_size,'fontname',Font_name);
% h_x=xlabel('$\bar{P}$','Interpreter','LaTex');
% set(h_x,'position',[3.3   -0.2    2.0001]);
% ylabel('$\gamma_R$','Interpreter','Latex');
% 
% 
% if print_fig~=0
%     
%     resolution='-r120';
%     
%     print(fig_mes_diag,'-zbuffer','-depsc',resolution,'m_es_Phase_Diagram.eps');
%     print(fig_mes_diag,'-zbuffer','-dpng',resolution,'m_es_Phase_Diagram.png');
%     
%     close(fig_mes_diag);
% end





















