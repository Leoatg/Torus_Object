

%m_ds vs m_es  R0 dependence

tic;
ua=7.7e-3;
GammaC=1;
R1=8.4e-3;

P1=1.9;
GammaR=1.4;

R0_min=43;
R0_space=0.2;
R0_max=45;

R0_series=R0_min:R0_space:R0_max;

print_fig=0;

Font_size=22;
Font_name='Times New Roman';

dk=0.1;
k=0:dk:100; %k axis
N_k=length(k);

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

m_min=0;
m_max=90;

N_R0=length(R0_series);
%N_GaR=length(GaR_series);

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

m_ds_R0=zeros(1,N_R0);
m_ds_R0(:,:)=m_max;


for nR0=1:N_R0
    
    R0=R0_series(nR0);
    
    m_ds_found_flag=0;
    
    g=2*ua;
    n0=Gammac/R1;
    Pth=(GammaR*Gammac)/R1;
    Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));
    
    for m=m_min:m_max
        
        parfor it=1:N_k
            
            L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1);
                -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1);
                -1i*GammaC*Phi0^2, -1i*GammaC*Phi0^2, -1i*(GammaR+R1*Phi0^2)];
            
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
        for it=2:length(Im_omega_t)
            if Im_omega_t(it)>0
                m_ds_R0(1,nR0)=m;
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


%fing plot range
N_R0=length(R0_series);

Nds=N_R0;
for itt=1:length(R0_series)
    if m_ds_R0(1,itt)==100
        Nds=itt;
        break;
    end
end

fig_m_ds_m_es_R0=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
%plot(R0_series(1:Nds),m_ds_R0(1,1:Nds),'linewidth',2);
h_mds=area(R0_series(1:Nds),m_ds_R0(1,1:Nds),'FaceColor',[0.5,0.9,0.6],'EdgeColor',[0,0.5,0.1]);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('R_0');
%ylabel('$m_{ds}$','interpreter','latex');
title(sprintf('GammaR=%.2f,  P1=%.2f',GammaR,P1));

hold all;
%plot(R0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds));
h_mes=area(R0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds),'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
hold off;

%h=legend('$m_{ds}$','$m_{es}$');
%set(h,'interpreter','latex');
%set(h,'position',[0.2554    0.8189    0.1312    0.0767]);

h=legend('Dynamically Stable','Energetically Satable');
set(h,'position',[0.2235    0.7973    0.3625    0.1200]);


if print_fig~=0
    
    print(fig_m_ds_m_es_R0,'-painter','-depsc','-r250',sprintf('m_ds_m_es_R0_Dependence_GammaR_P1=%.1f_P1=%.1f_Elena.eps',GammaR,P1));
    
    close(fig_m_ds_m_es_R0);
    
end











