

%m_ds vs m_es  P1 dependence

tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;


%GammaR=0.8;
%R0=28*sqrt(5/2);
R0=25;
P1=2.1;

GaR_min=0.8;
GaR_space=0.1;
GaR_max=2;

GaR_series=GaR_min:GaR_space:GaR_max;

print_fig=0;

Font_size=20;
Font_name='Times New Roman';

dk=0.1;
k=0:dk:80; %k axis
N_k=length(k);

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

m_min=0;
m_max=80;

N_GaR=length(GaR_series);
%N_GaR=length(GaR_series);

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

m_ds_GaR=zeros(1,N_GaR);
m_ds_GaR(:,:)=m_max;


for nGaR=1:N_GaR
    
    GammaR=GaR_series(nGaR);
    
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
        for it=2:length(Im_omega_t)
            if Im_omega_t(it)>0
                m_ds_GaR(1,nGaR)=m;
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
N_GaR=length(GaR_series);

Nds=N_GaR;
for itt=1:length(GaR_series)
    if m_ds_GaR(1,itt)==100
        Nds=itt;
        break;
    end
end

fig_m_ds_m_es_GaR=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
%plot(R0_series(1:Nds),m_ds_R0(1,1:Nds),'linewidth',2);
h_mds=area(GaR_series(1:Nds),m_ds_GaR(1,1:Nds),'FaceColor',[0.5,0.9,0.6],'EdgeColor',[0,0.5,0.1]);
axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
xlabel('GammaR');
%ylabel('$m_{ds}$','interpreter','latex');
ylabel('m');
title(sprintf('P1=%.2f,  R0=%.2f',P1,R0));

hold all;
%plot(R0_series(1:Nds),sqrt(ua*GammaR/R1)*R0_series(1:Nds));
h_mes=area(GaR_series(1:Nds),sqrt(ua*GaR_series*(P1-1)/R1)*R0,'FaceColor',[0.7,0.7,0.7],'EdgeColor','k');
hold off;

%h=legend('$m_{ds}$','$m_{es}$');
%set(h,'interpreter','latex');
%set(h,'position',[0.2554    0.8189    0.1312    0.0767]);

h=legend('Dynamically Stable','Energetically Satable');
set(h,'position',[0.4397    0.1640    0.3625    0.1200]);


if print_fig~=0
    
    print(fig_m_ds_m_es_GaR,'-painter','-depsc2','-r250',sprintf('m_ds_m_es_GammaR_Dependence_P1=%.1f_R0=%.1f.eps',P1,R0));
    print(fig_m_ds_m_es_GaR,'-painter','-dpng','-r250',sprintf('m_ds_m_es_GammaR_Dependence_P1=%.1f_R0=%.1f.png',P1,R0));
    
    %close(fig_m_ds_m_es_P1);
    
end











