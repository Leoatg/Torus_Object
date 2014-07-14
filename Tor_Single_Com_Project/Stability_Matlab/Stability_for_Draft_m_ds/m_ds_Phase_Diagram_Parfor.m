
% m_es phase diagram
tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

%R0=16*sqrt(5/2);
R0=28*sqrt(5/2);

P1_min=1.1;
P1_space=0.1;
P1_max=5+P1_space;

P1_series=P1_min:P1_space:P1_max;


GaR_min=0.1;
GaR_space=0.1;
GaR_max=4+GaR_space;

GaR_series=GaR_min:GaR_space:GaR_max;
 
if length(P1_series)~=length(GaR_series)
    fprintf('Error: mesh lengh not fit\n');
    fprintf('N P1=%d, N_GammaR=%d\n',length(P1_series),length(GaR_series));
    return;
end

print_fig=1;

Font_size=22;
Font_name='Times New Roman';

dk=0.01;
k=0:dk:100; %k axis
N_k=length(k);

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

m_min=0;
m_max=110;

N_P1=length(P1_series);

m_ds=zeros(N_P1,N_P1);
m_ds(:,:)=m_max;

n_P1=0;
n_Ga=0;

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle


parfor n_P1=1:N_P1
    
    P1=P1_series(n_P1);
    
    for n_Ga=1:N_P1
        GammaR=GaR_series(n_Ga);
        
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
%                 Eigen_omega_2(it)=complex(omega_tp(2,1),omega_tp(2,2));
%                 Eigen_omega_3(it)=complex(omega_tp(3,1),omega_tp(3,2));
            end
            
            %finding critical m by positive Imag omega
            Im_omega_t=imag(Eigen_omega_1);
            for it=floor(10/dk):length(Im_omega_t)
                if Im_omega_t(it)>0
                    break;
                end
            end
            if it~=length(Im_omega_t) || ((it==length(Im_omega_t)) && Im_omega_t(length(Im_omega_t))>0)
                m_ds(n_Ga,n_P1)=m;
                fprintf(sprintf('P1=%.2f, GammaR=%.2f, mds=%d\n',P1,GammaR,m));
                m_ds_found_flag=1;
                break;
            end
 
        end
        
        if m_ds_found_flag~=0
            continue;
        end
        
    end
    
end



[x,y]=meshgrid(P1_series,GaR_series);


fig_mds_diag=figure('position',[300 200 1000 800],'renderer','zbuffer','paperpositionmode','auto');

mesh(x,y,m_ds/R0);
axis square;
set(gca,'xlim',[P1_min,P1_max],'ylim',[GaR_min,GaR_max]);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
view(0,90);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
h_x=xlabel('$\bar{P}$','Interpreter','LaTex');
%set(h_x,'position',[2.1963   -0.0093    2.0001],'fontsize',16);
ylabel('$\gamma_R$','Interpreter','Latex','fontsize',16);
zlabel('$m_{ds}$','Interpreter','LaTex');
%title('$m_{ds}/R_0$','Interpreter','Latex');

%interp
[X,Y]=meshgrid(P1_min:0.001:P1_max,GaR_min:0.001:GaR_max);
[m_ds_inted]=interp2(x,y,m_ds,X,Y);

fig_mds_diag_interped=figure('renderer','zbuffer','paperpositionmode','auto');
mesh(X,Y,m_ds_inted/R0);
axis square;
set(gca,'xlim',[P1_min,P1_max],'ylim',[GaR_min,GaR_max]);
set(gca,'fontsize',Font_size,'Fontname',Font_name);
view(0,90);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
h_x=xlabel('$\bar{P}$','Interpreter','LaTex');
%set(h_x,'position',[2.1963   -0.0093    2.0001],'fontsize',16);
ylabel('$\gamma_R$','Interpreter','Latex','fontsize',16);
%title('$m_{ds}/R_0$','Interpreter','Latex');

if print_fig~=0
    
    resolution='-r150';
    
    % print(fig_mes_diag,'-zbuffer','-depsc',resolution,'m_ds_Phase_Diagram.eps');
    print(fig_mds_diag,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram.png');
    print(fig_mds_diag_interped,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram.png');
    print(fig_mds_diag_interped,'-zbuffer','-depsc',resolution,'m_ds_Phase_Diagram.eps');
    
    %close(fig_mes_diag);
end

save('Latest_Data.mat');

fprintf('Program finished.\n');

toc;



















