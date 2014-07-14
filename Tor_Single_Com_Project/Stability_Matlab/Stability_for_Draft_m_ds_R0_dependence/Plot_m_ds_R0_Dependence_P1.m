
%Plot_m_ds_R0_Dependence_P1

print_fig=0;

Font_size=20;
Font_name='Times New Roman';

%fing plot range
N_R0=length(R0_series);

Nds=N_R0;
for itt=1:length(R0_series)
    if m_ds_R0(1,itt)==100
        Nds=itt;
        break;
    end
end

fig_m_ds_R0=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
plot(R0_series(1:Nds),m_ds_R0(1,1:Nds),'linewidth',2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('R_0');
ylabel('$m_{ds}$','interpreter','latex');
%title(gca,sprintf('Gamma_R=%.1f, P bar=%.1f',GammaR,P1));
hold all;
st_legend=['$\bar{P}=',sprintf('%.1f',P1_series(1)),'$'];
for it=2:N_P1
    Nds=N_R0;
    for itt=1:length(R0_series)
        if m_ds_R0(it,itt)==100
            Nds=itt;
            break;
        end
    end
    
    plot(R0_series(1:Nds),m_ds_R0(it,1:Nds),'linewidth',2);
    st_legend=char(st_legend,['$\bar{P}=',sprintf('%.1f',P1_series(it)),'$']);
end
hold off;
h_l_P1=legend(st_legend,'Location','NorthEastOutside');
set(h_l_P1,'interpreter','latex');
set(h_l_P1,'fontsize',16);
set(h_l_P1,'position',[0.7 0.15 0.08 0.25]);


if print_fig~=0
    
    print(fig_m_ds_R0,'-painter','-depsc','-r250',sprintf('m_ds_R0_Dependence_P1_GammaR=%.1f.eps',GammaR));
    
    close(fig_m_ds_R0);
    
end