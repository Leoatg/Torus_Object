

fig_m_ds_m_es_R0=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
%plot(R0_series(1:Nds),m_ds_R0(1,1:Nds),'linewidth',2);
h_mds=area(R0_series(1:Nds)/sqrt(5/2),m_ds_R0(1,1:Nds),'FaceColor',[0.5,0.9,0.6],'EdgeColor',[0,0.5,0.1]);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('w_0');
%ylabel('$m_{ds}$','interpreter','latex');
title(sprintf('GammaR=%.2f,  P1=%.2f',GammaR,P1));

hold all;
%plot(R0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds));
h_mes=area(R0_series(1:Nds)/sqrt(5/2),sqrt(ua*GammaR/R1)*R0_series(1:Nds),'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
hold off;

%h=legend('$m_{ds}$','$m_{es}$');
%set(h,'interpreter','latex');
%set(h,'position',[0.2554    0.8189    0.1312    0.0767]);

h=legend('Dynamically Stable','Energetically Satable');
set(h,'position',[0.2235    0.7973    0.3625    0.1200]);