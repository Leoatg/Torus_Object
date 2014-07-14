
% m_es phase diagram
tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

R0=16*sqrt(5/2);

P1_min=1.1;
P1_space=0.2;
P1_max=3+P1_space;

P1_series=P1_min:P1_space:P1_max;


GaR_min=0.1;
GaR_space=0.2;
GaR_max=2+GaR_space;

GaR_series=GaR_min:GaR_space:GaR_max;
 
if length(P1_series)~=length(GaR_series)
    fprintf('Error: mesh lengh not fit\n');
    fprintf('N P1=%d, N_GammaR=%d\n',length(P1_series),length(GaR_series));
    return;
end

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

N_P1=length(P1_series);

m_ds=zeros(N_P1,N_P1);
m_ds(:,:)=m_max;

n_P1=0;
n_Ga=0;

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle


for P1=P1_min:P1_space:P1_max
    
    n_P1=n_P1+1;
    
    n_Ga=0;
    for GammaR=GaR_min:GaR_space:GaR_max
        
        n_Ga=n_Ga+1;
        
        m_ds_found_flag=0;

        g=2*ua;
        n0=Gammac/R1;
        Pth=(GammaR*Gammac)/R1;
        Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));
        
        for m=m_min:m_max
            
            parfor it=1:N_k
                
                L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
                    -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
                    -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
                
                
                omega_temp=eig(L,'nobalance');
%                 omega_de(:,1)=real(omega_temp);
%                 omega_de(:,2)=imag(omega_temp);
%                 omega_de(:,3)=abs(omega_temp);
%                 omega_de(:,4)=angle(omega_temp);
                omega_tp=zeros(3,2);
                omega_tp(:,1)=real(omega_temp);
                omega_tp(:,2)=imag(omega_temp);
                omega_tp=sortrows(omega_tp,-2);
                %if m~=0
                    %omega_de=sortrows(omega_de,-2);%sortrow 4 times
                %end
%                 Eigen_omega_1(it)=complex(omega_de(1,1),omega_de(1,2));
%                 Eigen_omega_2(it)=complex(omega_de(2,1),omega_de(2,2));
%                 Eigen_omega_3(it)=complex(omega_de(3,1),omega_de(3,2));
                  Eigen_omega_1(it)=complex(omega_tp(1,1),omega_tp(1,2));
                  Eigen_omega_2(it)=complex(omega_tp(2,1),omega_tp(2,2));
                  Eigen_omega_3(it)=complex(omega_tp(3,1),omega_tp(3,2));
            end
            
            %finding critical m by positive Imag omega
            Im_omega_t=imag(Eigen_omega_1);
            for itt=floor(10/dk):length(Im_omega_t)
                if Im_omega_t(itt)>0
                    break;
                end
            end
            if itt~=length(Im_omega_t) || ((itt==length(Im_omega_t)) && Im_omega_t(length(Im_omega_t))>0)
                m_ds(n_P1,n_Ga)=m;
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



[x,y]=meshgrid(GaR_series,P1_series);


fig_mds_diag=figure('position',[300 200 1000 800],'renderer','zbuffer','paperpositionmode','auto');

mesh(x,y,m_ds);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
view(0,90);
colormap(jet(4096));
shading interp;
colorbar('fontsize',Font_size,'fontname',Font_name);
h_x=ylabel('$\bar{P}$','Interpreter','LaTex');
%set(h_x,'position',[2.1963   -0.0093    2.0001],'fontsize',16);
xlabel('$\gamma_R$','Interpreter','Latex','fontsize',16);



if print_fig~=0
    
    resolution='-r150';
    
    % print(fig_mes_diag,'-zbuffer','-depsc',resolution,'m_es_Phase_Diagram.eps');
    print(fig_mes_diag,'-zbuffer','-dpng',resolution,'m_es_Phase_Diagram.png');
    
    %close(fig_mes_diag);
end


fprintf('Program finished.\n');

toc;



















