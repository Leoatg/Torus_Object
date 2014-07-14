

% Stability for Draft: R_0 dependent

tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

%R0=16*sqrt(5/2);

%P1=2;
GammaR=1.2;

R0_min=15;
R0_space=0.1;
R0_max=35;

R0_series=R0_min:R0_space:R0_max;

P1_min=1.8;
P1_space=0.2;
P1_max=3;

P1_series=P1_min:P1_space:P1_max;

% GaR_min=1;
% GaR_space=0.1;
% GaR_max=1.3;
% 
% GaR_series=GaR_min:GaR_space:GaR_max;

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
m_max=100;

N_R0=length(R0_series);
N_P1=length(P1_series);

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

m_ds_R0=zeros(N_P1,N_R0);
m_ds_R0(:,:)=m_max;

for n_P1=1:N_P1
    
    P1=P1_series(n_P1);
    
    parfor nR0=1:N_R0
        
        R0=R0_series(nR0);
        
        m_ds_found_flag=0;
        Eigen_omega_1=complex(zeros(1,N_k));
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
                    m_ds_R0(n_P1,nR0)=m;
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
st_legend=sprintf('P1=%.2f',P1_series(1));
for it=2:N_P1
    plot(R0_series,m_ds_R0(it,:),'linewidth',2);
    st_legend=char(st_legend,sprintf('P1=%.2f',P1_series(it)));
end
hold off;
legend(st_legend,'Location','NorthEastOutside');


if print_fig~=0
    
    resolution='-r150';
    
    %print(fig_m_ds_R0,'-painter','-depsc',resolution,'m_ds_R0_Dependence_P1.eps');
    print(fig_m_ds_R0,'-painter','-dpng','-r250','m_ds_R0_Dependence_P1.png');
    
    
    close(fig_m_ds_R0);
    
end

save(sprintf('m_ds_R0_Dependence_P1_GammaR=%.1f.mat',GammaR));

toc;

















