
%Solve k_star

global m;

Font_size=22;
Font_name='Times New Roman';

print_fig=0;

m_series=60:1:80;
N_m=length(m_series);
%m=60;

k_star1=zeros(1,N_m);
k_star2=k_star1;

P=2.5;
GammaR=1.4;
w0=16;
R0=w0*sqrt(5/2);

%k*1
for m_n=1:N_m
    
    m=m_series(m_n);
    
    options = optimset('Display','iter','FunValCheck','on');
    
    k1 = fzero(@f_Re1,50,options);
    
    k2=  fzero(@f_Re1,70,options);
    
    k_star1(1,m_n)=k1;
    k_star2(1,m_n)=k2;
    
    
end

N_k2_end=N_m;
for it=N_m:-1:1
    if (k_star2(it)-k_star1(it))>1
    N_k2_end=it;
    break;
    end
end

fig_k_star=figure('position',[300 200 1000 800],'renderer','painter','paperpositionmode','auto');
plot(m_series,k_star1,'linewidth',2);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
hold all;
plot(m_series(1:N_k2_end),k_star2(1:N_k2_end),'linewidth',2);
hold off;
h_l=legend('$k^*_1$','$k^*_2$');
set(h_l,'interpreter','latex','position',[0.6410 0.6125 0.1299 0.2585],'fontsize',Font_size);


if print_fig~=0
    
  print(fig_k_star,'-painter','-depsc','-r300',['k^star_',sprintf('P1=%.2f,GammaR=%.2f,R0=%.1f',P,GammaR,R0),'_Solved.eps']);
  close(fig_k_star);
end




