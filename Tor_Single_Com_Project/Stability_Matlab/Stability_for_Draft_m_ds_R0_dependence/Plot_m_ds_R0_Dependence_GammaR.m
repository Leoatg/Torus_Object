
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
st_legend=['$\gamma_R=',sprintf('%.1f',GaR_series(1)),'$'];
for it=2:N_GaR
    Nds=N_R0;
    for itt=1:N_R0
        if m_ds_R0(it,itt)==100
            Nds=itt;
            break;
        end
    end
    
    plot(R0_series(1:Nds),m_ds_R0(it,1:Nds),'linewidth',2);
    st_legend=char(st_legend,['$\gamma_R=',sprintf('%.1f',GaR_series(it)),'$']);
end
hold off;
h_l_Gr=legend(st_legend,'Location','NorthEastOutside');
set(h_l_Gr,'interpreter','latex');
set(h_l_Gr,'fontsize',16);
set(h_l_Gr,'position',[0.7 0.165 0.08 0.25]);


if print_fig~=0
    
    print(fig_m_ds_R0,'-painter','-depsc','-r250',sprintf('m_ds_R0_Dependence_GammaR_P1=%.1f.eps',P1));
    
    close(fig_m_ds_R0);
    
end





