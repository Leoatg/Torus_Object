
%plot k_star overdamped

load('P0=2.5,GammaR=1.4.mat');


print_fig=1;

Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;

%plot preparation
k_star=squeeze(k_star_all(1,1,:,:));
Nk1=1;
for it=1:length(k_star(2,:))
    if k_star(1,it)==0
    Nk1=it-1;
    break;
    end
end
Nk2_a=1;
Nk2_b=length(k_star(2,:));
for it=1:length(k_star(2,:))
    if k_star(2,it)~=0
    Nk2_a=it;
    break;
    end
end
for it=length(k_star(2,:)):-1:1
    if k_star(2,it)~=0
    Nk2_b=it;
    break;
    end
end


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

fig_k_star_interped=figure('position',[300 200 windows_width windows_hight],'renderer','painter','paperpositionmode','auto');
plot(k_star(3,1:Nk1),k_star(1,1:Nk1),'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
set(gca,'xlim',[32 90],'ylim',[10 100]);
xlabel('m');
hold all;
plot(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),'linewidth',3,'linestyle','--');
plot(m_axis_1(1:N_exp1),k1(1:N_exp1),'linewidth',2.5,'linestyle','-.','Color','r');
plot(m_axis_2(1:N_exp2),k2(1:N_exp2),'linewidth',2.5,'linestyle','-.','Color','r');
hold off;
% h_l=legend('$k^*_1$','$k^*_2$');
% set(h_l,'interpreter','latex','position',[0.6410 0.6125 0.1299 0.2585],'fontsize',Font_size);

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(c)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.82    0.08    0.08]);

if print_fig~=0
    
    
    print(fig_k_star_interped,'-painter','-depsc2','-r300',['k^star_',sprintf('P1=%.2f,GammaR=%.2f_Overdamped',P1_max,GaR_max),'.eps']);
  
    
    close(fig_k_star_interped);
end




