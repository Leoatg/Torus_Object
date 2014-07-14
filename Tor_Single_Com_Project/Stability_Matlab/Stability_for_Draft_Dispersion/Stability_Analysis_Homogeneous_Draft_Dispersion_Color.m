

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
dk=0.1;
k=-90:dk:90; %k axis

k_max=max(k);
k_min=min(k);

print_fig=0;

Font_size=24;
Font_name='Times New Roman';

Line_width=3;

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
        
        [V,D]=eig(L,'nobalance');
        
        
        v1=V(1:3,1).';
        omega_temp(1,1:3)=v1; omega_temp(1,4)=D(1,1);omega_temp(1,5)=real(v1(1,3));
        
        v2=V(1:3,2).';
        omega_temp(2,1:3)=v2; omega_temp(2,4)=D(2,2);omega_temp(2,5)=real(v2(1,3));
        
        v3=V(1:3,3).';
        omega_temp(3,1:3)=v3; omega_temp(3,4)=D(3,3);omega_temp(3,5)=real(v3(1,3));
        
        omega_temp=sortrows(omega_temp,5);
        
        omega_sorted_1(it)=omega_temp(1,4);
        omega_sorted_2(it)=omega_temp(2,4);
        omega_sorted_3(it)=omega_temp(3,4);
        
    end
    
    real_1=real(omega_sorted_1);imag_1=imag(omega_sorted_1);
    real_2=real(omega_sorted_2);imag_2=imag(omega_sorted_2);
    real_3=real(omega_sorted_3);imag_3=imag(omega_sorted_3);
    
    %----------------------------------------------------------
    %     %sort omega by hand
    if it_P_series==1 %Underdamped
        %real part
        %swap points
        [real_3(1),real_1(1)]=swap_num(real_3(1),real_1(1));
        [real_3(2),real_2(1)]=swap_num(real_3(2),real_2(1));
        %swap 1-2
        curve_temp=real_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)<0
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=real_1(swap_n1:swap_n2);
        real_1(swap_n1:swap_n2)=real_2(swap_n1:swap_n2);
        real_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        %swap 1-3
        curve_temp=real_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_1(swap_n1:end);
        real_1(swap_n1:end)=real_3(swap_n1:end);
        real_3(swap_n1:end)=temp_data;
        clear temp_data
        
        %single point adjust
        for it=floor(N_k/2)+1:N_k-1 %adjust 1
            if real_1(it+1)-real_1(it)>0.1
               real_1(it+1)=real_1(it);
               break;
            end
        end

        for it=10:-1:1 %adjust 3
            if real_3(it+1)-real_3(it)>0.1
               real_3(it)=real_3(it+1);
            end
        end
        %----------------------------------------------------
        %imag part
        imag_1(1)=imag_3(1);
        imag_2(1)=imag_3(2);
        imag_3(1:2)=imag_3(3);
        %swap 1-2
        curve_temp=imag_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.1
                if i_n<N_k/2
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=imag_1(swap_n1:swap_n2);
        imag_1(swap_n1:swap_n2)=imag_2(swap_n1:swap_n2);
        imag_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        %swap 1-3
        swap_n1=floor(N_k/2+1);
        temp_data=imag_1(swap_n1:end);
        imag_1(swap_n1:end)=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=temp_data;
        clear temp_data
        
    elseif it_P_series==2 %Bov
        %real part
        %swap 1-2
        curve_temp=real_2;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)<0
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=real_1(swap_n1:swap_n2);
        real_1(swap_n1:swap_n2)=real_2(swap_n1:swap_n2);
        real_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        
        %swap points
        [real_3(floor(N_k/2)+1),real_1(floor(N_k/2)+1)]=swap_num(real_3(floor(N_k/2)+1),real_1(floor(N_k/2)+1));
        [real_3(floor(N_k/2)-1),real_1(floor(N_k/2)-1)]=swap_num(real_3(floor(N_k/2)-1),real_1(floor(N_k/2)-1));
        
        %single point adjust
        for it=floor(N_k/2)+1:N_k-1 %adjust 3
            if real_3(it+1)-real_3(it)>0.1
               real_3(it+1)=real_3(it);
               break;
            end
        end
        
         for it=20:-1:1 %adjust 3
            if real_3(it+1)-real_3(it)>0.1
               real_3(it)=real_3(it+1);
            end
         end
        
        for it=10:-1:1 %adjust 1
            if real_1(it+1)-real_1(it)>0.1
               real_1(it)=real_1(it+1);
            end
        end
        for it=10:-1:1 %adjust 2
            if real_2(it+1)-real_2(it)>0.1
               real_2(it)=real_2(it+1);
            end
        end
         
        
        %imag part
        %swap 1-2
        curve_temp=imag_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.1
                if i_n<N_k/2
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=imag_1(swap_n1:swap_n2);
        imag_1(swap_n1:swap_n2)=imag_2(swap_n1:swap_n2);
        imag_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        
        %swap 1-3
        swap_n1=floor(N_k/2+1);
        temp_data=imag_1(swap_n1:end);
        imag_1(swap_n1:end)=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=temp_data;
        clear temp_data
        
    elseif it_P_series==3 %overdamped
        
        %real part
        %swap 1-2
        curve_temp=real_2;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)<0
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=real_1(swap_n1:swap_n2);
        real_1(swap_n1:swap_n2)=real_2(swap_n1:swap_n2);
        real_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        
        %single point adjust
        for it=floor(N_k/2)+1:N_k-1
            if real_3(it+1)-real_3(it)>0.1
               real_3(it+1)=real_3(it);
               break;
            end
        end
        
         for it=20:-1:1
            if real_3(it+1)-real_3(it)>0.1
               real_3(it)=real_3(it+1);
            end
         end
        
        for it=20:-1:1 %adjust 3
            if real_3(it+1)-real_3(it)>0.1
               real_3(it)=real_3(it+1);
            end
         end
        
        for it=10:-1:1 %adjust 1
            if real_1(it+1)-real_1(it)>0.1
               real_1(it)=real_1(it+1);
            end
        end
        for it=10:-1:1 %adjust 2
            if real_2(it+1)-real_2(it)>0.1
               real_2(it)=real_2(it+1);
            end
        end
        
         
         
        
        %imag part
        curve_temp=imag_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.1
                if i_n<N_k/2
                    swap_n1=i_n+1;
                else
                    swap_n2=i_n;
                end
            end
        end
        temp_data=imag_1(swap_n1:swap_n2);
        imag_1(swap_n1:swap_n2)=imag_2(swap_n1:swap_n2);
        imag_2(swap_n1:swap_n2)=temp_data;
        clear temp_data
        
        
    end
    
    %---------------------------------------------------------
    %plot Re_omega 
    
    c_value=['b','g','m']; %set color
    
    plot_k_displace=0.1;
    
    fig_Re=figure('position',[50 200 800 600],'renderer','painters','paperpositionmode','auto');
    plot(k,real_1,'LineWidth',Line_width,'color',c_value(1));
    axis_re=gca;
    if it_P_series==1
        ax_po=get(gca,'position');
    end
    set(gca,'ActivePositionProperty','position');
    %axis square;
    set(gca,'fontsize',Font_size,'fontname',Font_name);
    set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    set(gca,'xlim',[k_min+plot_k_displace,k_max-plot_k_displace]);
    set(gca,'ylim',[-5,5]);
    if it_P_series==3
        xlabel('k','fontsize',30);
    end
    hold(axis_re,'all');
    plot(axis_re,k,real_2,'LineWidth',Line_width,'color',c_value(2));
    plot(axis_re,k,real_3,'LineWidth',Line_width,'color',c_value(3));
    plot(axis_re,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_re,'off');
    
    
    %plot Im_omega
    
    %adjust Im limit
    Im_max=max([max(imag_1),max(imag_2),max(imag_3)]);
    if Im_max<0.5
        Im_max=0.5;
    end
    Im_min=min([min(imag_1),min(imag_2),min(imag_3)]);
    exceed_ratio_down=1.1;
    exceed_ratio_up=1;
    
    c_value=['m','g','b']; %set color
    
    fig_Im=figure('position',[800 200 800 600],'renderer','painters','paperpositionmode','auto');
    plot(gca,k,imag_1,'LineWidth',Line_width,'color',c_value(1));
    axis_im=gca;
    set(gca,'ActivePositionProperty','position');
    %axis square;
    set(gca,'fontsize',Font_size,'fontname',Font_name);
    set(gca,'xlim',[k_min+plot_k_displace,k_max-plot_k_displace]);
    set(gca,'ylim',[Im_min*exceed_ratio_down,Im_max*exceed_ratio_up]);
    set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    y_tick=get(gca,'ytick');
    if y_tick<0.5
    y_tick=[y_tick, 0.5];
    end
    set(gca,'ytick',y_tick);
    set(gca,'yticklabel',[]);
    set(gca,'yticklabel',num2str(sortrows(y_tick.',1),'%.1f'));
    if it_P_series==3
        xlabel('k','fontsize',30);
    end
    ylabel('Im \omega');
    hold(axis_im,'all');
    plot(axis_im,k,imag_2,'LineWidth',Line_width,'color',c_value(2));
    plot(axis_im,k,imag_3,'LineWidth',Line_width,'color',c_value(3));
    plot(axis_im,k,zeros(1,N_k),'linestyle','-.','linewidth',2,'color','r');
    hold(axis_im,'off');
    
    
    %plot subplot
    m=Parameter_series(it_P_series+3,1);
    
    for it=1:N_k
        
        L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
            -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
            -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
        
        [V,D]=eig(L,'nobalance');
        
        
        v1=V(1:3,1).';
        omega_temp(1,1:3)=v1; omega_temp(1,4)=D(1,1);omega_temp(1,5)=real(v1(1,3));
        
        v2=V(1:3,2).';
        omega_temp(2,1:3)=v2; omega_temp(2,4)=D(2,2);omega_temp(2,5)=real(v2(1,3));
        
        v3=V(1:3,3).';
        omega_temp(3,1:3)=v3; omega_temp(3,4)=D(3,3);omega_temp(3,5)=real(v3(1,3));
        
        omega_temp=sortrows(omega_temp,1);
        
        omega_sorted_1(it)=omega_temp(1,4);
        omega_sorted_2(it)=omega_temp(2,4);
        omega_sorted_3(it)=omega_temp(3,4);
        
    end
    
    real_1=real(omega_sorted_1);imag_1=imag(omega_sorted_1);
    real_2=real(omega_sorted_2);imag_2=imag(omega_sorted_2);
    real_3=real(omega_sorted_3);imag_3=imag(omega_sorted_3);
    
    
    %--------------------------------------------------------------------------
    %     %sort omega by hand
    if it_P_series==1 %Underdamped
        %real part
        %swap 1-2
        curve_temp=real_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_1(swap_n1:end);
        real_1(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        %swap 3-2
        curve_temp=real_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_3(swap_n1:end);
        real_3(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        
        %imag part
        %swap 3-2
        curve_temp=imag_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=imag_2(swap_n1:end);
        imag_2(swap_n1:end)=temp_data;
        clear temp_data
        
        %    %swap 1-3
        curve_temp=imag_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.1
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=imag_1(swap_n1:end);
        imag_1(swap_n1:end)=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=temp_data;
        clear temp_data
        
    elseif it_P_series==2
        %real part
        %swap 1-2
        curve_temp=real_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_1(swap_n1:end);
        real_1(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        %swap 3-2
        curve_temp=real_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_3(swap_n1:end);
        real_3(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        
        %imag part
        %swap 3-2
        curve_temp=imag_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=imag_2(swap_n1:end);
        imag_2(swap_n1:end)=temp_data;
        clear temp_data
        
        %    %swap 1-3
        curve_temp=imag_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.1
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=imag_1(swap_n1:end);
        imag_1(swap_n1:end)=imag_3(swap_n1:end);
        imag_3(swap_n1:end)=temp_data;
        clear temp_data
        
    elseif it_P_series==3 %overdamped
        %real part
        %swap 1-2
        curve_temp=real_1;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_1(swap_n1:end);
        real_1(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        %swap 3-2
        curve_temp=real_3;
        for i_n=1:N_k-1
            if abs(curve_temp(i_n+1)-curve_temp(i_n))>0.5
                if curve_temp(i_n+1)
                    swap_n1=i_n+1;
                end
            end
        end
        temp_data=real_3(swap_n1:end);
        real_3(swap_n1:end)=real_2(swap_n1:end);
        real_2(swap_n1:end)=temp_data;
        clear temp_data
        
        
        
        
    end
    %---------------------------------------------------------------------------
    %Subplot Re_omega
    
    if it_P_series==1%set line color
        c_value=['b','m','g'];
    elseif it_P_series==2
        c_value=['b','m','g'];
    elseif it_P_series==3
        c_value=['b','m','g'];
    end
    
    figure(fig_Re);
    axis_sub_fig_Re=axes('Position',[0.6100    0.2100    0.2800    0.2800]);
    plot(axis_sub_fig_Re,k,real_1,'LineWidth',Line_width,'color',c_value(1));
    set(axis_sub_fig_Re,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    set(axis_sub_fig_Re,'xlim',[0,k_max]);
    set(axis_sub_fig_Re,'ylim',[-5,5]);
    %xlabel(axis_sub_fig_Re,'k');
    %ylabel(axis_sub_fig_Re,'Re \omega');
    %title(axis_re,sprintf('m=%d, Gamma_R=%.1f, P bar=%.1f',m,GammaR,P1));
    hold(axis_sub_fig_Re,'all');
    plot(axis_sub_fig_Re,k,real_2,'LineWidth',Line_width,'color',c_value(2));
    plot(axis_sub_fig_Re,k,real_3,'LineWidth',Line_width-1,'color',c_value(3));
    % plot(axis_sub_fig_Re,k,zeros(1,N_k),'linestyle','-.','linewidth',1,'color','r');
    hold(axis_sub_fig_Re,'off');
    
    
    %Subplot Im_omega
    
    %adjust Im limit
    Im_max=max([max(imag_1),max(imag_2),max(imag_3)]);
    if Im_max<0.5
        Im_max=0.5;
    end
    Im_min=min([min(imag_1),min(imag_2),min(imag_3)]);
    exceed_ratio_down=1.1;
    exceed_ratio_up=1;
       
    c_value=['b','m','g']; %set color
  
    
    figure(fig_Im);
    axis_sub_fig_Im=axes('Position',[0.6100    0.2100    0.2800    0.2800]);
    plot(axis_sub_fig_Im,k,imag_1,'LineWidth',Line_width,'color',c_value(1));
    set(axis_sub_fig_Im,'xlim',[0,k_max]);
    set(gca,'ylim',[Im_min*exceed_ratio_down,Im_max*exceed_ratio_up]);
    set(axis_sub_fig_Im,'PlotBoxAspectRatio',PlotBoxAspectRatio);
    y_tick=get(gca,'ytick');
    if y_tick<0.5
    y_tick=[y_tick, 0.5];
    end
    set(gca,'ytick',y_tick);
    set(gca,'yticklabel',[]);
    set(gca,'yticklabel',num2str(sortrows(y_tick.',1),'%.1f'));
    %xlabel(axis_sub_fig_Im,'k');
    % ylabel(axis_sub_fig_Im,'Im \omega');
    hold(axis_sub_fig_Im,'all');
    plot(axis_sub_fig_Im,k,imag_3,'LineWidth',Line_width,'color',c_value(3));
    plot(axis_sub_fig_Im,k,imag_2,'LineWidth',Line_width-1,'color',c_value(2));
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





