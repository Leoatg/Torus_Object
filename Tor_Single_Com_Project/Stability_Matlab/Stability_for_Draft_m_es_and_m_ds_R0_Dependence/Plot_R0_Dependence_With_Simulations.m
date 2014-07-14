
%Plot R0 Dependance with simulations

ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

P1=2.107;
GammaR=1.4;

w0_min=20;
w0_space=0.2;
w0_max=30;

w0_series=w0_min:w0_space:w0_max;


R0_series=w0_series*sqrt(5/2);

print_fig=1;

Font_size=22;
Font_name='Times New Roman';


Dy_sim_w0_m=[20,58;
             22,63;
             24,68;
             26,74;
             28,79];
         
Ey_sim_w0_m=[20,44;
             22,49;
             24,54;
             26,60;
             28,65];
         
Stable_sim_w0_m=[20,25;
                 22,29;
                 24,33;
                 26,36;
                 28,40];


%plot using R0
fig_m_ds_m_es_R0_sim=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
h_mds=area(R0_series(1:Nds),m_ds_R0(1,1:Nds),'FaceColor',[0.5,0.9,0.6],'EdgeColor',[0,0.5,0.1]);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('$r_0$','interpreter','latex');
%title(sprintf('GammaR=%.2f,  P1=%.2f',GammaR,P1));
hold all;
h_mes=area(R0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds),'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
marker_size=120;
scatter(Dy_sim_w0_m(:,1)*sqrt(5/2),Dy_sim_w0_m(:,2),marker_size,'k','s','fill');
scatter(Ey_sim_w0_m(:,1)*sqrt(5/2),Ey_sim_w0_m(:,2),marker_size,'b','d','fill');
scatter(Stable_sim_w0_m(:,1)*sqrt(5/2),Stable_sim_w0_m(:,2),marker_size,'r','o','fill');
hold off;

h=legend('Dynamically stable','Energetically satable','Quick decay','Prolonged decay','Stable');
set(h,'position',[0.5100    0.1347    0.3725    0.3133]);

if print_fig~=0
    
    print(fig_m_ds_m_es_R0_sim,'-painter','-depsc','-r350',sprintf('m_ds_m_es_R0_Sim_Dependence_GammaR=%.1f_P1=%.3f.eps',GammaR,P1));
    
    %print(fig_m_ds_m_es_w0,'-painter','-depsc','-r250',sprintf('m_ds_m_es_w0_Dependence_GammaR=%.1f_P1=%.3f.eps',GammaR,P1));
    
    close(fig_m_ds_m_es_R0_sim);
    
end




%plot using w0
% fig_m_ds_m_es_R0=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
% h_mds=area(w0_series(1:Nds),m_ds_R0(1,1:Nds),'FaceColor',[0.5,0.9,0.6],'EdgeColor',[0,0.5,0.1]);
% axis square;
% set(gca,'fontsize',Font_size,'fontname',Font_name);
% xlabel('w_0');
% title(sprintf('GammaR=%.2f,  P1=%.2f',GammaR,P1));
% hold all;
% h_mes=area(w0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds),'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
% hold off;
% 
% h=legend('Dynamically Stable','Energetically Satable');
% set(h,'position',[0.2235    0.7973    0.3625    0.1200]);







