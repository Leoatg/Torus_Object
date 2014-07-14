

% stability analysis, dispertion analysis
% plot omega
GammaR=1.5;
P1=1.62;

ua=7.7e-3;
Gammac=1;
R1=8.4e-3;

g=2*ua;
n0=Gammac/R1;
Pth=(GammaR*Gammac)/R1;
Phi0=sqrt((1/R1)*(P1*Pth/n0-GammaR));

R0=16*sqrt(5/2);

m=0;%persistent current m
dk=0.01;

k=-100:dk:50; %k axis

print_fig=0;
Record_Mov_m=0;
m_max=70;

Fontsize=18;
N=length(k);

omega=(2*k*m*R1+sqrt(k.^4 * R1^2+4*k.^2*R0^2*R1*ua*GammaR*(P1-1)))/(2*R0^2*R1);
omega2=(2*k*m*R1-sqrt(k.^4 * R1^2+4*k.^2*R0^2*R1*ua*GammaR*(P1-1)))/(2*R0^2*R1);

Negative_area=0;
Negative_1=1e7;
for it=1:N
    if omega(it)<0
        Negative_1=it;
        break;
    end
end
if Negative_1<=N
    for it=Negative_1:N
        if omega(it)>0
            Negative_2=it;
            break;
        end
    end
    
    Negative_area=sum(omega(Negative_1:Negative_2))*dk;
end

fig_omega=figure('renderer','painter','paperpositionmode','auto');
plot(k,omega,'linewidth',2);
set(gca,'fontsize',Fontsize);
xlabel('k');
ylabel('\omega');
hold all
plot(k,zeros(1,N),'linewidth',2,'linestyle','--');
%plot(k,omega2,'linewidth',2);
hold off
st_title=sprintf('GammaR=%.1f P1=%.2f m=%d Narea=%.1f', GammaR, P1,m,Negative_area);
title(st_title);

%search negative point
for it=N:-1:1
    if omega(it)<0
        break;
    end
end
if it~=1
text(k(it),omega(it),'\bullet','fontsize',20,'HorizontalAlignment','center');
end

%print(gcf,'-painter','-dpng','-r250',sprintf('Dispersion curve_Overdamped_m=%d.png',m));



fprintf('\nProgram finished.\n');





