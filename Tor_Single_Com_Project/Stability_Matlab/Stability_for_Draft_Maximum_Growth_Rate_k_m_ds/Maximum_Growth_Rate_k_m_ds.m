
% dynamic instability maximum growth rate


ua=7.7e-3;
R1=8.4e-3;
GammaC=1;
g=2*ua;

P1=2.2;
GammaR=1.2;

R0=28*sqrt(5/2);

n0=GammaC/R1;
Phi0=sqrt(GammaR*(P1-1)/R1);

m_min=0;
m_max=120;
m_series=m_min:m_max;

Font_size=20;
Font_name='Times New Roman';

dk=0.1;
k=0:dk:100; %k axis
N_k=length(k);

k_begin=10;

print_fig=0;

N_m=length(m_series);

Imomega_min=zeros(1,N_m);
k_min_d=Imomega_min;

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

n_m=0;
for m=m_min:m_max
    
    n_m=n_m+1;
    
    for it=1:N_k
        
        L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
            -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
            -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
        
        omega_temp=eig(L,'nobalance');
        omega_tp=zeros(3,2);
        omega_tp(:,1)=real(omega_temp);
        omega_tp(:,2)=imag(omega_temp);
        omega_tp=sortrows(omega_tp,-2);
        %end
        Eigen_omega_1(it)=complex(omega_tp(1,1),omega_tp(1,2));
        
    end
    
    %finding critical m by positive Imag omega
    Im_omega_t=imag(Eigen_omega_1);
        
    [o_m,k_m]=max(Im_omega_t(floor(k_begin/dk):end));
    
    if o_m>0
       Imomega_min(n_m)=o_m;
       k_min_d(n_m)=k(k_m);
    end
    
    fprintf(sprintf('m=%d  Im omega=%.1f   k_d=%.1f\n',m,o_m,k(k_m)));
    
end

fig_o_min=figure('renderer','painter','paperpositionmode','auto');
plot(m_series,Imomega_min,'linewidth',2);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
axis square;
xlabel('m');
ylabel('min \omega');
st_title=sprintf('GammaR=%.1f, P1=%.1f, R_0=%.1f',GammaR,P1,R0);
title(st_title);

fig_k_min=figure('renderer','painter','paperpositionmode','auto');
plot(m_series,k_min_d,'linewidth',2);
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



















