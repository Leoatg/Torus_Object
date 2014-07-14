

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

m=55;%persistent current m
dk=0.001;

k=0:dk:80; %k axis

print_fig=0;

Font_size=22;
Font_name='Times New Roman';

N=length(k);
Eigen_omega_1_Re_sorted=complex(zeros(1,N));%three eigen omega(k)
Eigen_omega_2_Re_sorted=complex(zeros(1,N));
Eigen_omega_3_Re_sorted=complex(zeros(1,N));

Eigen_omega_1_Im_sorted=Eigen_omega_1_Re_sorted;
Eigen_omega_2_Im_sorted=Eigen_omega_1_Re_sorted;
Eigen_omega_3_Im_sorted=Eigen_omega_1_Re_sorted;

omega_temp=complex(zeros(3,1),0);
omega_de=zeros(3,4);%set by column: 1 real part; 2 imaginary part; 3 amplitude; 4 angle

omega_Re_sorted=omega_de;
omega_Im_sorted=omega_de;


Parameter_series=[0,0.5,2;
                   0 1.5 2.5;
                   60 1.5 2.5];
               
Parameter_series(:,2)=Parameter_series(:,2)*R0^2;
%Parameter_series(:,3)=Parameter_series(:,3)*R0^2;
               
               
PlotBoxAspectRatio=[1.6 1 1];

for it_P_series=1:3
   
    m=Parameter_series(it_P_series,1);
    GammaR=Parameter_series(it_P_series,2);
    P1=Parameter_series(it_P_series,3);

    ua=7.7e-3*R0^2;
    Gammac=1*R0^2;
    R1=8.4e-3*R0^2;
    
    g=2*ua;
    n0=Gammac/R1;
    Pth=(GammaR*Gammac)/(R1);
    Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));

for it=1:N
    
    L=[((m+k(it))^2-m^2)/(2)+ua*Phi0^2, ua*Phi0^2, (g+0.5i*R1)*Phi0;
        -ua*Phi0^2, -(((m-k(it))^2-m^2)/(2)+ua*Phi0^2), -(g-0.5i*R1)*Phi0;
        -1i*R1*Phi0*n0, -1i*R1*Phi0*n0, -1i*(GammaR+R1*Phi0^2)];
    
    
    omega_temp=eig(L,'nobalance')/R0^2;
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


fig_Re=figure('position',[50 200 800 600],'renderer','painter','paperpositionmode','auto');
plot(k,real(Eigen_omega_1_Re_sorted),'linewidth',3,'color','b');
axis_re=gca;
%axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
xlabel('k');
ylabel('Re \omega');
%title(axis_re,sprintf('m=%d, Gamma_R=%.1f, P bar=%.1f',m,GammaR,P1));
hold(axis_re,'all');
plot(axis_re,k,real(Eigen_omega_2_Re_sorted),'linewidth',3,'color','b');
plot(axis_re,k,real(Eigen_omega_3_Re_sorted),'linewidth',3,'color','b');
plot(axis_re,k,zeros(1,N),'linestyle','-.','linewidth',2,'color','r');
hold(axis_re,'off');

%Im_max=max(imag(Eigen_omega_1_Im_sorted));
Im_max=0.5;
Im_min=min([min(imag(Eigen_omega_1_Im_sorted)),min(imag(Eigen_omega_2_Im_sorted)),min(imag(Eigen_omega_3_Im_sorted))]);
exceed_ratio=1.2;

fig_Im=figure('position',[800 200 800 600],'renderer','painter','paperpositionmode','auto');
plot(gca,k,imag(Eigen_omega_1_Im_sorted),'linewidth',3,'color','b');
axis_im=gca;
%axis square;
set(gca,'fontsize',Font_size,'fontname',Font_name);
%set(gca,'ylim',[Im_min*exceed_ratio,Im_max]);
set(gca,'PlotBoxAspectRatio',PlotBoxAspectRatio);
% y_tick=get(gca,'ytick');
% if max(y_tick)<0.5
%    y_tick=[y_tick 0.5];
%    set(gca,'ytick',y_tick);
% end
xlabel('k');
ylabel('Im \omega');
hold(axis_im,'all');
plot(axis_im,k,imag(Eigen_omega_2_Im_sorted),'linewidth',3,'color','b');
plot(axis_im,k,imag(Eigen_omega_3_Im_sorted),'linewidth',3,'color','b');
plot(axis_im,k,zeros(1,N),'linestyle','-.','linewidth',2,'color','r');
hold(axis_im,'off');


if print_fig~=0 
    
 print(fig_Re,'-painter','-depsc','-r250',sprintf('Dispersion_Curve_Re_m=%d_GammaR=%.1f_P1=%.1f_try.eps',m,GammaR,P1));
 print(fig_Im,'-painter','-depsc','-r250',sprintf('Dispersion_Curve_Im_m=%d_GammaR=%.1f_P1=%.1f_try.eps',m,GammaR,P1));


close(fig_Re);
close(fig_Im);
end

    
end


fprintf('\nProgram finished.\n');





