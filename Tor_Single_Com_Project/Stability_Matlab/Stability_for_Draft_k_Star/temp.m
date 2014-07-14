%interpt
m_start=50.2;
m_axis_1=m_start:0.01:110;
m_axis_2=m_start:0.01:64;
k1=interp1(k_star(3,1:Nk1),k_star(1,1:Nk1),m_axis_1,'PCHIP','extrap');
k2=interp1(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),m_axis_2,'PCHIP','extrap');

N_exp1=1;
N_exp2=1;
for it=1:length(k1)
    if k1(it)<=k_star(1,1);
        N_exp1=it;
        break;
    end
end
for it=1:length(k1)
    if k2(it)>=k_star(2,Nk2_a);
        N_exp2=it;
        break;
    end
end

fig_k_star_interped=figure('position',[300 200 1000 800],'renderer','painter','paperpositionmode','auto');
plot(k_star(3,1:Nk1),k_star(1,1:Nk1),'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
set(gca,'xlim',[48,110]);
xlabel('m');
hold all;
plot(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),'linewidth',3,'linestyle','--');
plot(m_axis_1(1:N_exp1),k1(1:N_exp1),'linewidth',2.5,'linestyle',':','Color','r');
plot(m_axis_2(1:N_exp2),k2(1:N_exp2),'linewidth',2.5,'linestyle',':','Color','r');
hold off;
h_l=legend('$k^*_1$','$k^*_2$','Extrapolated');
set(h_l,'interpreter','latex','position',[0.6410 0.6125 0.1299 0.2585],'fontsize',Font_size);