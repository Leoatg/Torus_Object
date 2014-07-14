
% energy instability maximum growth rate


ua=7.7e-3;
R1=8.4e-3;
GammaC=1;

P1=2.5;
GammaR=1.5;

R0=16*sqrt(5/2);

n0=GammaC/R1;
Phi0=sqrt(GammaR*(P1-1)/R1);

k=-250:0.1:100;

m_series=0:100;

Font_size=20;
Font_name='Times New Roman';

print_fig=0;

N_m=length(m_series);

omega_min=zeros(1,N_m);
k_min=omega_min;

for it=40:N_m
    
    m=m_series(it);
    
    omega=(2*k*m+sqrt(k.^4+4*ua*R0^2*Phi0^2*k.^2))/(2*R0^2);
    
    [o_m,k_m]=min(omega);
    
    if o_m<0
       omega_min(it)=o_m;
       k_min(it)=k(k_m);
    end
end

fig_o_min=figure('renderer','painter','paperpositionmode','auto');
plot(m_series,omega_min,'linewidth',2);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
axis square;
xlabel('m');
ylabel('min \omega');
st_title=sprintf('GammaR=%.1f, P1=%.1f, R_0=%.1f',GammaR,P1,R0);
title(st_title);

fig_k_min=figure('renderer','painter','paperpositionmode','auto');
plot(m_series,k_min,'linewidth',2);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
ylabel('k min');
title(st_title);


% if print_fig~=0
%     
%     print(fig_o_min,'-zbuffer','-depsc','-r250','m_es_Phase_Diagram.eps');
%     print(fig_k_min,'-zbuffer','-dpng','-r250','m_es_Phase_Diagram.png');
%     
%     close(fig_mes_diag);
% end



















