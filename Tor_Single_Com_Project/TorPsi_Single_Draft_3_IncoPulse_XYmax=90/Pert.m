
clear fre_series

theta=atan2(y,x);

Pert_Strength=0.005;

phi0=10;

u_0=1;
v_0=1;
w_0=1;


k_per=50;

Pert_psi=u_0*exp(1i*(phi0+k_per)*theta)+conj(v_0)*exp(1i*(phi0-k_per)*theta);
%Pert_psi=Pert_psi+Noise;

%Pert_n=(w_0*exp(1i*k_per*theta)+w_0*exp(-1i*k_per*theta)+Noise);

k_min=2;
k_space=5;
k_max=40;

N_k1=length(k_min:k_space:k_max);
k_phase_shift=2*pi/N_k1;

% figure('renderer','zbuffer');
% mesh(x,y,abs(Pert_psi));
% axis square;
% view(0,90);
% colorbar;
% title('pert');

for k1=k_min:k_space:k_max
    
    k_phase_shift=2*pi/k1;
    
    Pert_psi=Pert_psi+exp(1i*(phi0+k1+k_phase_shift)*theta)+exp(1i*(phi0-k1-k_phase_shift)*theta);
    %Pert_n=Pert_n+exp(1i*(k1+k_phase_shift)*theta)+exp(1i*(k1-k_phase_shift)*theta);
    
%     mesh(x,y,abs(Pert_psi));
%     axis square;
%     view(0,90);
%     colorbar;
%     title('pert');
%     
%     a=1;
    
end

psi_mode=abs(psi);

Pert_psi=psi_mode.*Pert_psi*Pert_Strength;

%fre_series=cos((m0+k)*theta)+cos((m0-k)*theta);

%fre_series=sqrt((cos(phi_1)+cos(phi_2)).^2+(sin(phi_2)-sin(phi_1)).^2);

% for it=10:40:50
%
% end

%pert=psi_mode.*fre_series*0.1;


% figure('renderer','zbuffer');
% mesh(x,y,abs(psi_mode+Pert_psi).^2);
% axis square;
% view(0,90);
% colorbar;


figure('renderer','zbuffer','position',[100 200 1600 800],'paperpositionmode','auto');

subplot(1,2,1);
mesh(x,y,abs(Pert_psi));
set(gca,'fontsize',20);
axis square;
set(gca,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
%colorbar;
xlabel('x');
ylabel('y');
title('perturbation');

subplot(1,2,2)
mesh(x,y,abs(Pert_psi));
set(gca,'fontsize',20);
axis square;
set(gca,'xlim',[-XYmax,XYmax],'ylim',[-XYmax,XYmax]);
view(0,90);
colorbar;
xlabel('x');
ylabel('y');


print(gcf,'-zbuffer','-dpng','Perturbation.png');




