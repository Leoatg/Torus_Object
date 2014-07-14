

% stability analysis, dispertion analysis

% P1=2.5;
% GammaR=1.5;

ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

% g=2*ua;
% n0=Gammac/R1;
% Pth=(GammaR*Gammac)/R1;
% Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));

R0=16*sqrt(5/2);
%R0=30;

%dk=0.001;
dk=0.01;
k=-100:dk:100; %k axis

k_max=max(k);
k_min=min(k);

print_fig=0;

Font_size=24;
Font_name='Times New Roman';

N_k=length(k);
Eigen_omega_1_Re_sorted=complex(zeros(1,N_k));%three eigen omega(k)
Eigen_omega_2_Re_sorted=complex(zeros(1,N_k));
Eigen_omega_3_Re_sorted=complex(zeros(1,N_k));

Eigen_omega_1_Im_sorted=Eigen_omega_1_Re_sorted;
Eigen_omega_2_Im_sorted=Eigen_omega_1_Re_sorted;
Eigen_omega_3_Im_sorted=Eigen_omega_1_Re_sorted;

omega_temp=complex(zeros(5,3),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

omega_sorted_1=complex(zeros(1,N_k));
omega_sorted_2=complex(zeros(1,N_k));
omega_sorted_3=complex(zeros(1,N_k));

% omega_Re_sorted=omega_de;
% omega_Im_sorted=omega_de;

m=60;%persistent current m
m_mid=30;
P0=2.5;
Ga_bov=4*(P0-1)/P0^2;

Parameter_series=[m   0.5   P0;        %m \gammaR P1
                  m Ga_bov  P0;
                  m   1.4   P0; 
                  0   0.5   P0;
                  0 Ga_bov  P0;
                  0   1.4   P0];
              


PlotBoxAspectRatio=[1.6 1 1];

for it_P_series=1:3
    
    
    GammaR=Parameter_series(it_P_series,2);
    P1=Parameter_series(it_P_series,3);
    
    g=2*ua;
    n0=Gammac/R1;
    Pth=(GammaR*Gammac)/R1;
    Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));
    
  
    
    m=Parameter_series(it_P_series,1);
    
    for it=1:N_k
        
        L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
            -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
            -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
        
        %[V,D]=eig(L,'nobalance');
        [V,D]=eig(L);
        
        v1=V(1:3,1).'; Mo1=sqrt(abs(v1));
        omega_temp(1,1:3)=v1; omega_temp(1,4)=D(1,1);omega_temp(1,5)=real(v1(1,3))/Mo1;
        
        v2=V(1:3,2).';Mo2=real(dot(v2,conj(v2)));
        omega_temp(2,1:3)=v2; omega_temp(2,4)=D(2,2);omega_temp(2,5)=real(v2(1,3))/Mo2;
        
        v3=V(1:3,3).';Mo3=real(dot(v3,conj(v3)));
        omega_temp(3,1:3)=v3; omega_temp(3,4)=D(3,3);omega_temp(3,5)=real(v3(1,3))/Mo3;
        
        omega_temp=sortrows(omega_temp,5);
         
        omega_sorted_1(it)=Mo1;%omega_temp(1,1);
        omega_sorted_2(it)=Mo2;%omega_temp(2,1);
        omega_sorted_3(it)=Mo3;%omega_temp(3,1);
        
    end
    
    %plot Re_omega
    fig_Re=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
    plot(k,real(omega_sorted_1),'linewidth',3,'color','b');
    axis_re=gca;
    if it_P_series==1
        ax_po=get(gca,'position');
    end
   set(gca,'ActivePositionProperty','position');
    %axis square;
    set(gca,'fontsize',Font_size,'fontname',Font_name);
    set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    set(gca,'xlim',[k_min,k_max]);
    if it_P_series==3
    xlabel('k','fontsize',30);
    end
    %ylabel('Re \omega');
    %title(axis_re,sprintf('m=%d, Gamma_R=%.1f, P bar=%.1f',m,GammaR,P1));
    hold(axis_re,'all');
    plot(axis_re,k,real(omega_sorted_2),'linewidth',3,'color','g');
    plot(axis_re,k,real(omega_sorted_3),'linewidth',3,'color','m');
    plot(axis_re,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_re,'off');

   
    %Im_max=max(imag(Eigen_omega_1_Im_sorted));
%     Im_max=0.4;
%     Im_min=min([min(imag(Eigen_omega_1_Im_sorted)),min(imag(Eigen_omega_2_Im_sorted)),min(imag(Eigen_omega_3_Im_sorted))]);
%     exceed_ratio=1.2;
    
    %plot Im_omega
    fig_Im=figure('position',[800 200 800 600],'renderer','painter','paperpositionmode','auto');
    plot(gca,k,imag(omega_sorted_1),'linewidth',3,'color','b');
    axis_im=gca;
  set(gca,'ActivePositionProperty','position');
    %axis square;
    set(gca,'fontsize',Font_size,'fontname',Font_name);
    set(gca,'xlim',[k_min,k_max]);
    %set(gca,'ylim',[Im_min*exceed_ratio,Im_max]);
    set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    y_tick=get(gca,'ytick');
    if max(y_tick)<0.4
        y_tick=[y_tick 0.4];
        set(gca,'ytick',y_tick);
    end
    y_tick=get(gca,'ytick');
    set(gca,'yticklabel',[]);
    set(gca,'yticklabel',num2str(sortrows(y_tick.'),'%.1f'));
    if it_P_series==3
    xlabel('k','fontsize',30);
    end
    %ylabel('Im \omega');
    hold(axis_im,'all');
    plot(axis_im,k,imag(omega_sorted_2),'linewidth',3,'color','g');
    plot(axis_im,k,imag(omega_sorted_3),'linewidth',3,'color','m');
    plot(axis_im,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_im,'off');
    
    %plot subplot
    m=Parameter_series(it_P_series+3,1);
    
    for it=1:N_k
        
        L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
            -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
            -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
        
        
        omega_temp=eig(L,'nobalance');
        omega_de(:,1)=real(omega_temp);
        omega_de(:,2)=imag(omega_temp);
        omega_de(:,3)=abs(omega_temp);
        omega_de(:,4)=angle(omega_temp);
        
        
        omega_Re_sorted=sortrows(omega_de,-1);
        
        omega_Im_sorted=sortrows(omega_de,-2);
        
        
        Eigen_omega_1_Re_sorted(it)=complex(omega_Re_sorted(1,1),omega_Re_sorted(1,2));
        Eigen_omega_2_Re_sorted(it)=complex(omega_Re_sorted(2,1),omega_Re_sorted(2,2));
        Eigen_omega_3_Re_sorted(it)=complex(omega_Re_sorted(3,1),omega_Re_sorted(3,2));
        
        Eigen_omega_1_Im_sorted(it)=complex(omega_Im_sorted(1,1),omega_Im_sorted(1,2));
        Eigen_omega_2_Im_sorted(it)=complex(omega_Im_sorted(2,1),omega_Im_sorted(2,2));
        Eigen_omega_3_Im_sorted(it)=complex(omega_Im_sorted(3,1),omega_Im_sorted(3,2));
        
    end
    
    
    %--------------------------------------------------------------------------
    
    %Subplot Re_omega
    figure(fig_Re);
    axis_sub_fig_Re=axes('Position',[0.6100    0.2100    0.2800    0.2800]);
    plot(axis_sub_fig_Re,k,real(Eigen_omega_1_Re_sorted),'linewidth',3,'color','b');
    %set(axis_sub_fig_Re,'fontsize',Font_size,'fontname',Font_name);
    set(axis_sub_fig_Re,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    set(axis_sub_fig_Re,'xlim',[0,k_max]);
    %xlabel(axis_sub_fig_Re,'k');
    %ylabel(axis_sub_fig_Re,'Re \omega');
    %title(axis_re,sprintf('m=%d, Gamma_R=%.1f, P bar=%.1f',m,GammaR,P1));
    hold(axis_sub_fig_Re,'all');
    plot(axis_sub_fig_Re,k,real(Eigen_omega_2_Re_sorted),'linewidth',3,'color','b');
    plot(axis_sub_fig_Re,k,real(Eigen_omega_3_Re_sorted),'linewidth',3,'color','b');
    plot(axis_sub_fig_Re,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_sub_fig_Re,'off');
    
    
    %plot Im_omega
    figure(fig_Im);
    axis_sub_fig_Im=axes('Position',[0.6100    0.2100    0.2800    0.2800]);
    plot(axis_sub_fig_Im,k,imag(Eigen_omega_1_Im_sorted),'linewidth',3,'color','b');
    %set(gca,'fontsize',Font_size,'fontname',Font_name);
    set(axis_sub_fig_Im,'xlim',[0,k_max]);
    %set(axis_sub_fig_Im,'ylim',[Im_min*exceed_ratio,Im_max]);
    set(axis_sub_fig_Im,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    y_tick=get(axis_sub_fig_Im,'ytick');
%     if max(y_tick)<0.5
%         y_tick=[y_tick 0.5];
%         set(axis_sub_fig_Im,'ytick',y_tick);
%     end
    y_tick=get(axis_sub_fig_Im,'ytick');
    set(axis_sub_fig_Im,'yticklabel',[]);
    set(axis_sub_fig_Im,'yticklabel',num2str(sortrows(y_tick.'),'%g'));
     %xlabel(axis_sub_fig_Im,'k');
    % ylabel(axis_sub_fig_Im,'Im \omega');
    hold(axis_sub_fig_Im,'all');
    plot(axis_sub_fig_Im,k,imag(Eigen_omega_2_Im_sorted),'linewidth',3,'color','b');
    plot(axis_sub_fig_Im,k,imag(Eigen_omega_3_Im_sorted),'linewidth',3,'color','b');
    plot(axis_sub_fig_Im,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_sub_fig_Im,'off');
    
    
    
    
    
    
    
    
    %---------------------------------------------------------------------------
    if print_fig~=0
        
        currentFolder = pwd;
        
        if abs(GammaR-Ga_bov)<1e-2;
            dataFolder_name='Draft_Dispersion_Bov';
        elseif GammaR-Ga_bov>0
            dataFolder_name='Draft_Dispersion_Overdamped';
        elseif GammaR-Ga_bov<0
            dataFolder_name='Draft_Dispersion_Underdamped';
        end
        
        if exist(dataFolder_name,'dir')==0
            mkdir(dataFolder_name);
        end
       
        print(fig_Re,'-painters','-depsc2','-r350',[currentFolder,'\',dataFolder_name,'\',sprintf('Dispersion_Curve_Re_m=%d_GammaR=%.2f_P1=%.1f.eps',m,GammaR,P1)]);
        print(fig_Im,'-painters','-depsc2','-r350',[currentFolder,'\',dataFolder_name,'\',sprintf('Dispersion_Curve_Im_m=%d_GammaR=%.2f_P1=%.1f.eps',m,GammaR,P1)]);
        
        
        close(fig_Re);
        close(fig_Im);
        
    end
    
    
end


fprintf('\nProgram finished.\n');





