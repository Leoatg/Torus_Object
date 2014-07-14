

load('P0=2.5,GammaR=0.5.mat');

Font_size=22;
Font_name='Times New Roman';
windows_width=500;
windows_hight=500;

print_fig=1;

%plot preparation
k_star=squeeze(k_star_all(1,1,:,:));
%Nk1=1;
% for it=1:length(k_star(2,:))
%     if k_star(1,it)==0
%     Nk1=it-1;
%     break;
%     end
% end
Nk2_a=1;
for it=1:length(k_star(2,:))
    if k_star(2,it)>=1
       Nk2_a=it;
        break;
    end
end
Nk2_b=length(k_star(2,:));
% for it=1:length(k_star(2,:))
%     if k_star(2,it)~=0
%     Nk2_a=it;
%     break;
%     end
% end
for it=length(k_star(2,:)):-1:1
    if k_star(2,it)~=0
    Nk2_b=it;
    break;
    end
end

%%single plot
fig_k_star=figure('position',[300 200 windows_width windows_hight],'renderer','painter','paperpositionmode','auto');
plot([0,k_star(3,Nk2_a:Nk2_b)],[0.0001,k_star(1,Nk2_a:Nk2_b)],'linewidth',3);
axis square;
set(gca,'xlim',[0 k_star(3,Nk2_b)]);%,'ylim',[k_star(1,Nk2_a),k_star(1,Nk2_b)]);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
hold all;
plot([0,k_star(3,Nk2_a:Nk2_b)],[32.3,k_star(2,Nk2_a:Nk2_b)],'linewidth',3,'linestyle','--');
hold off;

h_l=legend('$k^*_1$','$k^*_2$');
set(h_l,'interpreter','latex','position',[0.6780    0.2    0.1299    0.2585],'fontsize',Font_size,'box','off');

h_a=annotation('textbox', [0.1460    0.750    0.08    0.08],'String','(a)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+2,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',[0.1460    0.82    0.08    0.08]);


if print_fig~=0

print(fig_k_star,'-painter','-depsc2','-r300',['k^star_',sprintf('P1=%.2f,GammaR=%.2f',P1_max,GaR_max),'_Underdamped.eps']);

close(fig_k_star);
end





