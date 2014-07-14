

% stability analysis, dispertion analysis

P1=2.1;
GammaR=1.4;%0.96;

ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

g=2*ua;
n0=Gammac/R1;
Pth=(GammaR*Gammac)/R1;
Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));

w0=22;
R0=w0*sqrt(5/2);

%m=80;%persistent current m
dk=0.1;

k=-350:dk:350; %k axis

print_fig=0;
Record_Mov_m=0;
show_fig=0;


Font_size=22;
Font_name='Times New Roman';

N=length(k);
Eigen_omega_1=complex(zeros(1,N));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N));
Eigen_omega_3=complex(zeros(1,N));

%energy second order perturbation for three modes
E2_1=zeros(1,N);
E2_2=E2_1;
E2_3=E2_1;

m_min=0;
m_max=80;

m_series=m_min:m_max;
N_m=length(m_series);

Imomega_max=zeros(1,N_m);%store maximum Im \omega
ImE2_min=Imomega_max;

% if Record_Mov_m~=0
%     writerObj = VideoWriter('Dispersion_mov');
%     writerObj.FrameRate=2;
%     open(writerObj);
%
%     m_min=50;
% else
%     m_min=m;
%     m_max=m;
% end

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

E2_temp=omega_temp;
E2_de=omega_de;

for m_it=1:N_m
    
    m=m_series(m_it);
    
    for it=1:N
        
        L=[((m+k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
            -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2*R0^2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
            -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
        
        omega_temp=eig(L,'nobalance');
        omega_de(:,1)=real(omega_temp);
        omega_de(:,2)=imag(omega_temp);
        omega_de(:,3)=abs(omega_temp);
        omega_de(:,4)=angle(omega_temp);
        if m~=0
            omega_de=sortrows(omega_de,-2);%sortrow for omega
        end
        Eigen_omega_1(it)=complex(omega_de(1,1),omega_de(1,2));
        Eigen_omega_2(it)=complex(omega_de(2,1),omega_de(2,2));
        Eigen_omega_3(it)=complex(omega_de(3,1),omega_de(3,2));
        
        [V,D] = eig(L,'nobalance');
        
        E2_1(it)=D(1,1)*(abs(V(1,1))^2-abs(V(2,1))^2)*2*pi+(real(V(1,1))+real(V(2,1)))*V(3,1)*pi*(g+0.5i*R1);
        E2_2(it)=D(2,2)*(abs(V(1,2))^2-abs(V(2,2))^2)*2*pi+(real(V(1,2))+real(V(2,2)))*V(3,2)*pi*(g+0.5i*R1);
        E2_3(it)=D(3,3)*(abs(V(1,3))^2-abs(V(2,3))^2)*2*pi+(real(V(1,3))+real(V(2,3)))*V(3,3)*pi*(g+0.5i*R1);
        
        E2_temp(1)=E2_1(it);
        E2_temp(2)=E2_2(it);
        E2_temp(3)=E2_3(it);
        E2_de(:,1)=real(E2_temp);
        E2_de(:,2)=imag(E2_temp);
        E2_de(:,3)=abs(E2_temp);
        E2_de(:,4)=angle(E2_temp);
        
        
            E2_de=sortrows(E2_de,-2);%sortrow for E2
        
        
        E2_1(it)=complex(E2_de(1,1),E2_de(1,2));
        E2_2(it)=complex(E2_de(2,1),E2_de(2,2));
        E2_3(it)=complex(E2_de(3,1),E2_de(3,2));
        
    end
    
    if m_it==1 && show_fig~=0
        %
        if m==m_min || Record_Mov_m==0
            fig_E2_im=figure('renderer','painter','position',[100 100 1800 800],'paperpositionmode','auto');
            %         fig_disp=figure('renderer','painters','position',[100 100 1600 800],'paperpositionmode','auto');
                     axis_re=subplot(1,2,1);
                     axis_im=subplot(1,2,2);
        end
        %
    end
    %
    %     axis_re=subplot(1,2,1);
    %     plot(axis_re,k,real(Eigen_omega_1),'linewidth',2);
    %     axis square;
    %     set(axis_re,'fontsize',Fontsize);
    %     xlabel('k');
    %     ylabel('Re \omega');
    %     title(axis_re,sprintf('m=%d, Gamma_R=%.2f, P bar=%.2f, R0=%.1f',m,GammaR,P1,R0));
    %     hold(axis_re,'all');
    %     plot(axis_re,k,real(Eigen_omega_2),'linewidth',2);
    %     plot(axis_re,k,real(Eigen_omega_3),'linewidth',2);
    %     plot(axis_re,k,zeros(1,N),'linestyle','-.','linewidth',2);
    %     hold(axis_re,'off');
    %     legend('\omega_1','\omega_2','\omega_3');
    %
    %     axis_im=subplot(1,2,2);
    %     plot(axis_im,k,imag(Eigen_omega_1),'linewidth',2);
    %     axis square;
    %     set(axis_im,'fontsize',Fontsize);
    %     xlabel('k');
    %     ylabel('Imag \omega');
    %     hold all
    %     plot(axis_im,k,zeros(1,N),'linestyle','-.','linewidth',2);
    %     plot(axis_im,k,imag(Eigen_omega_2),'linewidth',2);
    %     plot(axis_im,k,imag(Eigen_omega_3),'linewidth',2);
    %     hold off
    %     legend('\omega_1','\omega_2','\omega_3');
    if show_fig~=0
        
        subplot(1,2,1);
        plot(gca,k,real(E2_1),'linewidth',2);
        set(gca,'ylim',[-1,2]);
        set(gca,'fontsize',Font_size);
        xlabel('k');
        hold all
        plot(gca,k,real(E2_2),'linewidth',2);
        plot(gca,k,real(E2_3),'linewidth',2);
        hold off
        title('Re E2')
        h_l=legend('$E^{(2)} 1$','$E^{(2)} 2$','$E^{(2)} 3$');
        set(h_l,'Interpreter','latex');
        
        
        
        subplot(1,2,2);
        plot(gca,k,imag(E2_1),'linewidth',2);
        set(gca,'fontsize',Font_size);
        xlabel('k');
        hold all
        plot(gca,k,imag(E2_2),'linewidth',2);
        plot(gca,k,imag(E2_3),'linewidth',2);
        hold off
        title(['Im E2',sprintf('    m=%d',m)]);
        h_l=legend('$E^{(2)} 1$','$E^{(2)} 2$','$E^{(2)} 3$');
        set(h_l,'Interpreter','latex');
    
        pause(0.1);
    end
    
    
    %finding critical m by positive Imag omega
    Im_omega_t=imag(Eigen_omega_1);
    [omega_max,I]=max(Im_omega_t);
    %if omega_max>0
        Imomega_max(m_it)=omega_max;
        % h_t=text(k(I),omega_max,'\bullet','fontsize',20,'HorizontalAlignment','center');
    %end
    
    %find minimum E2
    Im_E2=imag(E2_3);
    [Im_E2_min,I]=min(Im_E2);
    ImE2_min(m_it)=Im_E2_min;
    % h_t=text(k(I),omega_max,'\bullet','fontsize',20,'HorizontalAlignment','center');
    
    if Record_Mov_m~=0
        writeVideo(writerObj,getframe(fig_disp));
    end
    
    
end

%plot Im omega and Im E2
fig_Imomega_max=figure('position',[300 400 800 600],'renderer','painter','paperpositionmode','auto');
plot(m_series,Imomega_max,'linewidth',3);
axis square;
set(gca,'fontsize',Font_size,'Fontname',Font_name);
xlabel('m');
%ylabel('max Im \omega');
hold all;
plot(m_series,abs(ImE2_min),'linewidth',3);
hold off;
h_l=legend('Im \omega max','Im E2 min (absolute value)');




%plot E^(2)
% fig_E2_re=figure('renderer','painter','position',[100 100+100 1600 800],'paperpositionmode','auto');
% plot(gca,k,real(E2_1),'linewidth',2);
% set(gca,'ylim',[-1,1]);
% set(gca,'fontsize',Fontsize);
% xlabel('k');
% hold all
% plot(gca,k,real(E2_2),'linewidth',2);
% plot(gca,k,real(E2_3),'linewidth',2);
% hold off
% title('Re E2')
% h_l=legend('$E^{(2)} 1$','$E^{(2)} 2$','$E^{(2)} 3$');
% set(h_l,'Interpreter','latex');

% figure;
% plot(gca,k,real(E2_1)+real(E2_2)+real(E2_3),'linewidth',2);
% title('Sum E2')

% fig_E2_im=figure('renderer','painter','position',[100 100 1600 800],'paperpositionmode','auto');
% plot(gca,k,imag(E2_1),'linewidth',2);
% set(gca,'fontsize',Fontsize);
% xlabel('k');
% hold all
% plot(gca,k,imag(E2_2),'linewidth',2);
% plot(gca,k,imag(E2_3),'linewidth',2);
% hold off
% title('Im E2')
% h_l=legend('$E^{(2)} 1$','$E^{(2)} 2$','$E^{(2)} 3$');
% set(h_l,'Interpreter','latex');


if print_fig~=0 && Record_Mov_m==0
    print(gcf,'-painters','-dpng','-r250',sprintf('E2_m=%d_GammaR=%.1f_P1=%.1f_w0=%d.png',m,GammaR,P1,w0));
    close(gcf);
end

if Record_Mov_m~=0
    close(writerObj);
end

fprintf('\nProgram finished.\n');





