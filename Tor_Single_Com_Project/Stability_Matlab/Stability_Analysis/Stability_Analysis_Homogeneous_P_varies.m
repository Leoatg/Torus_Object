

% stability analysis, dispertion analysis

%P1=2.1;
P1_series=2:0.01:2.8;
N_P1=length(P1_series);

GammaR=1.4;%0.96;

ua=7.7e-3;
Gammac=1;
R1=8.4e-3;
%R1=8.4e-10;

%g=2e-10*ua;
g=2*ua;
n0=Gammac/R1;
Pth=(GammaR*Gammac)/R1;
Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));

R0=24*sqrt(5/2);
%R0=57;

m=70;%persistent current m
dk=0.1;

k=-180:dk:180; %k axis

print_fig=0;
Record_Mov_m=1;
m_max=80;

Fontsize=18;


N=length(k);
Eigen_omega_1=complex(zeros(1,N));%three eigen omega(k)
Eigen_omega_2=complex(zeros(1,N));
Eigen_omega_3=complex(zeros(1,N));

if Record_Mov_m~=0
    writerObj = VideoWriter('Dispersion_P1_mov');
    writerObj.FrameRate=2;
    open(writerObj);
    
    m_min=50;
else
    m_min=m;
    m_max=m;
end

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle
for itn=1:N_P1
    
    P1=P1_series(itn);
    
    g=2*ua;
    n0=Gammac/R1;
    Pth=(GammaR*Gammac)/R1;
    Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));
    
    
    
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
    omega_de=sortrows(omega_de,-2);%sortrow 4 times
    end
    Eigen_omega_1(it)=complex(omega_de(1,1),omega_de(1,2));
    Eigen_omega_2(it)=complex(omega_de(2,1),omega_de(2,2));
    Eigen_omega_3(it)=complex(omega_de(3,1),omega_de(3,2));
end

if P1==min(P1_series) || Record_Mov_m==0
fig_disp=figure('renderer','painter','position',[100 100 1600 800],'paperpositionmode','auto');
axis_re=subplot(1,2,1);
axis_im=subplot(1,2,2);
end

axis_re=subplot(1,2,1);
plot(axis_re,k,real(Eigen_omega_1),'linewidth',2);
axis square;
set(axis_re,'fontsize',Fontsize);
xlabel('k');
ylabel('Re \omega');
title(axis_re,sprintf('m=%d, Gamma_R=%.2f, P bar=%.2f, R0=%.1f',m,GammaR,P1,R0));
hold(axis_re,'all');
plot(axis_re,k,real(Eigen_omega_2),'linewidth',2);
plot(axis_re,k,real(Eigen_omega_3),'linewidth',2);
plot(axis_re,k,zeros(1,N),'linestyle','-.','linewidth',2);
hold(axis_re,'off');


axis_im=subplot(1,2,2);
plot(axis_im,k,imag(Eigen_omega_1),'linewidth',2);
axis square;
set(axis_im,'fontsize',Fontsize,'ylim',[-5,0.5]);
xlabel('k');
ylabel('Imag \omega');
hold all
plot(axis_im,k,zeros(1,N),'linestyle','-.','linewidth',2);
plot(axis_im,k,imag(Eigen_omega_2),'linewidth',2);
plot(axis_im,k,imag(Eigen_omega_3),'linewidth',2);
hold off
%finding critical m by derivertive
% Im_omega_t=imag(Eigen_omega_2);
% Im_omega_t_diff=diff(Im_omega_t)/dk;
% for it=length(Im_omega_t_diff):-1:1
%     if abs(Im_omega_t_diff(it))<1e-5
%        break;
%     end
% end
% text(k(it),Im_omega_t(it),'\bullet','fontsize',20,'HorizontalAlignment','center');

%finding critical m by positive Imag omega
Im_omega_t=imag(Eigen_omega_1);
Im_omega_t_diff=diff(Im_omega_t)/dk;
for it=2:length(Im_omega_t)
    if Im_omega_t(it)>0
       break;
    end
end
if it~=length(Im_omega_t) || ((it==length(Im_omega_t)) && Im_omega_t(length(Im_omega_t))>0)
text(k(it),Im_omega_t(it),'\bullet','fontsize',20,'HorizontalAlignment','center');
title(axis_im,sprintf('k*1=%.2f',k(it)));
 end
% diff_k_axis=k(1:length(Im_omega_t_diff));
% figure;
% plot(diff_k_axis,Im_omega_t_diff);

if Record_Mov_m~=0
writeVideo(writerObj,getframe(fig_disp));
end

end

if print_fig~=0 && Record_Mov_m==0
print(gcf,'-painter','-dpng','-r250',sprintf('Dispersion_Curve_m=%d_GammaR=%.1f_P1=%.1f_R0=%.1f.png',m,GammaR,P1,R0));
close(fig_disp);
end

if Record_Mov_m~=0
    close(writerObj);
end

fprintf('\nProgram finished.\n');





