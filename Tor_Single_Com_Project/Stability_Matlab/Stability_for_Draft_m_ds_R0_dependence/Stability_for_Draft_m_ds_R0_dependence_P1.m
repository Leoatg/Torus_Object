

% Stability for Draft: R_0 dependent

tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

%R0=16*sqrt(5/2);

P1=2;
GammaR=1;

R0_min=23;
R0_space=0.1;
R0_max=35;

R0_series=R0_min:R0_space:R0_max;

GaR_min=1;
GaR_space=0.1;
GaR_max=1.3;

GaR_series=GaR_min:GaR_space:GaR_max;

print_fig=0;

Font_size=18;
Font_name='Times New Roman';

dk=0.1;
k=0:dk:20; %k axis
N_k=length(k);

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

m_min=0;
m_max=80;

N_R0=length(R0_series);
N_GaR=length(GaR_series);

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

m_ds_R0=zeros(N_GaR,N_R0);
m_ds_R0(:,:)=m_max;

for n_GaR=1:N_GaR
    
    GammaR=GaR_series(n_GaR);

for nR0=1:N_R0
    
    R0=R0_series(nR0);
    
    m_ds_found_flag=0;
    
    g=2*ua;
    n0=Gammac/R1;
    Pth=(GammaR*Gammac)/R1;
    Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));
    
    for m=m_min:m_max
        
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
        for it=floor(16/dk):length(Im_omega_t)
            if Im_omega_t(it)>0
                m_ds_R0(n_GaR,nR0)=m;
                fprintf(sprintf('R0=%.2f, P1=%.2f, GammaR=%.2f, mds=%d\n',R0,P1,GammaR,m));
                m_ds_found_flag=1;
                break;
            end
        end
        
        if m_ds_found_flag~=0
            break;
        end
        
    end
    
end

end

fig_m_ds_R0=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
plot(R0_series,m_ds_R0(1,:),'linewidth',2);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('R_0');
ylabel('$m_{ds}$','interpreter','latex');
%title(gca,sprintf('Gamma_R=%.1f, P bar=%.1f',GammaR,P1));
hold all;
st_legend=sprintf('GammaR=%.2f',GaR_series(1));
for it=2:N_GaR
    plot(R0_series,m_ds_R0(it,:),'linewidth',2);
    st_legend=char(st_legend,sprintf('GammaR=%.2f',GaR_series(it)));
end
hold off;
legend(st_legend,'Location','NorthEastOutside');


if print_fig~=0 
    
    resolution='-r150';
    
%print(fig_m_ds_R0,'-painter','-depsc',resolution,sprintf('Dispersion_Curve_Re_GammaR=%.1f_P1=%.1f.eps',GammaR,P1));
print(fig_m_ds_R0,'-painter','-dpng','-r250',sprintf('Dispersion_Curve_Im_GammaR=%.1f_P1=%.1f.png',GammaR,P1));


close(fig_m_ds_R0);

end























