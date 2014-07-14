
% m_es phase diagram
tic;
ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

R0=16*sqrt(5/2);

P1_min=2.5;
P1_space=0.1;
P1_max=2.5;

P1_series=P1_min:P1_space:P1_max;

P0=2.5;

Ga_bov=4*(P0-1)/P0^2;

GaR_min=1.4;
GaR_space=0.1;
GaR_max=GaR_min;

GaR_series=GaR_min:GaR_space:GaR_max;
%GaR_series=[0.5, Ga_bov, 1.5];

 
% if length(P1_series)~=length(GaR_series)
%     fprintf('Error: mesh lengh not fit\n');
%     fprintf('N P1=%d, N_GammaR=%d\n',length(P1_series),length(GaR_series));
%     return;
% end

print_fig=1;

Font_size=22;
Font_name='Times New Roman';

dk=0.001;
k=0:dk:100; %k axis
N_k=length(k);

Eigen_omega_1=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N_k));
Eigen_omega_3=complex(zeros(1,N_k));

m_min=0;
m_max=110;

N_P1=length(P1_series);
N_Ga=length(GaR_series);

m_ds=zeros(N_P1,N_P1);
m_ds(:,:)=m_max;

% n_P1=0;
% n_Ga=0;

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

N_m=length(m_min:m_max);

% k_star=zeros(3,2); %first three row k*1,k*2,m
% star_flag=0;
k_star_all=zeros(N_P1,N_Ga,3,N_m);
% k_star(1:2,:)=-1;
% k_star(3,:)=m_min:m_max;

for n_P1=1:N_P1
    
    P1=P1_series(n_P1);
    
    for n_Ga=1:N_Ga
        GammaR=GaR_series(n_Ga);
        k_star=zeros(3,N_m);
        star_flag=1;
        
        m_ds_found_flag=0;
        Eigen_omega_1=complex(zeros(1,N_k));
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
            for it=floor(5/dk):length(Im_omega_t)
                if Im_omega_t(it)>0
                    k_star(1,star_flag)=k(it);
                    k_star(3,star_flag)=m;
                    star_flag=star_flag+1;
                    fprintf('P1=%.1f, GammaR=%.1f, k1=%.1f,m=%d\n',P1,GammaR,k(it),m);
                    break;
                end
            end
            
            k1_it=it;
            
            for it=(k1_it+1):length(Im_omega_t)
                if Im_omega_t(it)<0
                    k_star(2,star_flag)=k(it);
                    break;
                end
            end
            
            k_star_all(n_P1,n_Ga,:,:)=k_star;
            
        end

    end
    
end

%plot preparation
k_star=squeeze(k_star_all(1,1,:,:));
Nk1=1;
for it=1:length(k_star(2,:))
    if k_star(1,it)==0
    Nk1=it-1;
    break;
    end
end
Nk2_a=1;
Nk2_b=length(k_star(2,:));
for it=1:length(k_star(2,:))
    if k_star(2,it)~=0
    Nk2_a=it;
    break;
    end
end
for it=length(k_star(2,:)):-1:1
    if k_star(2,it)~=0
    Nk2_b=it;
    break;
    end
end

%%single plot
fig_k_star=figure('position',[300 200 1000 800],'renderer','painter','paperpositionmode','auto');
plot(k_star(3,1:Nk1),k_star(1,1:Nk1),'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
hold all;
plot(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),'linewidth',3,'linestyle','--');
hold off;

%interpt
m_start=52.7;
m_axis_1=m_start:0.01:110;
m_axis_2=m_start:0.01:64;
k1=interp1(k_star(3,1:Nk1),k_star(1,1:Nk1),m_axis_1,'PCHIP','extrap');
k2=interp1(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),m_axis_2,'PCHIP','extrap');

N_exp1=1;
N_exp2=1;
for it=1:length(k1)
    if k1(it)<=k_star(1,1);
        N_exp1=it;
        break;
    end
end
for it=1:length(k1)
    if k2(it)>=k_star(2,Nk2_a);
        N_exp2=it;
        break;
    end
end

fig_k_star_interped=figure('position',[300 200 1000 800],'renderer','painter','paperpositionmode','auto');
plot(k_star(3,1:Nk1),k_star(1,1:Nk1),'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
hold all;
plot(k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b),'linewidth',3,'linestyle','--');
plot(m_axis_1(1:N_exp1),k1(1:N_exp1),'linewidth',2.5,'linestyle',':','Color','r');
plot(m_axis_2(1:N_exp2),k2(1:N_exp2),'linewidth',2.5,'linestyle',':','Color','r');
hold off;
h_l=legend('$k^*_1$','$k^*_2$','Extrapolated');
set(h_l,'interpreter','latex','position',[0.6410 0.6125 0.1299 0.2585],'fontsize',Font_size);

% plot(m_axis_1,k1,'linewidth',3);
% axis square;
% set(gca,'fontsize',Font_size,'Fontname',Font_name);
% xlabel('m');
% hold all;
% plot(m_axis_2,k2,'linewidth',3,'linestyle','--');
% hold off;



%%plotyy
% fig_k_star=figure('position',[300 200 1000 800],'renderer','painter','paperpositionmode','auto');
% [AX,H1,H2]=plotyy(k_star(3,1:Nk1),k_star(1,1:Nk1),k_star(3,Nk2_a:Nk2_b),k_star(2,Nk2_a:Nk2_b));
% axis square;
% set(AX(1),'fontsize',Font_size,'Fontname',Font_name);
% set(AX(2),'fontsize',Font_size,'Fontname',Font_name,'xtick',[]);
% xlabel(AX(1),'m');
% ylabel(AX(1),'$k_1^*$','interpreter','latex');
% ylabel(AX(2),'$k_2^*$','interpreter','latex');
% set(H1,'linewidth',2);
% set(H2,'linewidth',2);
% 
% title(AX(1),['$\bar{P}=$',sprintf('%.1f',P1_max),' $\gamma_R=$',sprintf('%.2f',GaR_max)],'interpreter','latex');

if print_fig~=0
    
    resolution='-r150';
    %print(fig_k_star,'-painter','-dpng','-r250',['k^star_',sprintf('P1=%.2f,GammaR=%.2f',P1_max,GaR_max),'.png']);
    %print(fig_k_star,'-painter','-depsc','-r250',['k^star_',sprintf('P1=%.2f,GammaR=%.2f',P1_max,GaR_max),'.eps']);
    print(fig_k_star_interped,'-painter','-depsc','-r300',['k^star_',sprintf('P1=%.2f,GammaR=%.2f',P1_max,GaR_max),'_interped.eps']);
    % print(fig_mes_diag,'-zbuffer','-depsc',resolution,'m_ds_Phase_Diagram.eps');
%     print(fig_mds_diag,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram.png');
%     print(fig_mds_diag_interped,'-zbuffer','-dpng',resolution,'m_ds_Phase_Diagram.png');
%     print(fig_mds_diag_interped,'-zbuffer','-depsc',resolution,'m_ds_Phase_Diagram.eps');
    
    %close(fig_mes_diag);
end

%save('Latest_Data.mat');

fprintf('Program finished.\n');

toc;



















