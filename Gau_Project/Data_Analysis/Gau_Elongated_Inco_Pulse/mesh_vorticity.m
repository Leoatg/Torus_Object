
%plot all

zoom_value=3;
Font_size=16;

figure('position',[100 100 1400 900],'paperpositionmode','auto','renderer','zbuffer');

subplot(2,3,1);
mesh(x,y,abs(psiL_temp).^2);
axis square;
set(gca,'fontsize',Font_size);
view(0,90);
colorbar;
zoom(zoom_value);
xlabel('x');
ylabel('y');
title('|\psi_L|^2');

subplot(2,3,2);
mesh(x,y,abs(psiR_temp).^2);
axis square;
set(gca,'fontsize',Font_size);
view(0,90);
colorbar;
zoom(zoom_value);
xlabel('x');
ylabel('y');
title('|\psi_R|^2');

subplot(2,3,3);
plot(Time,Energy);
axis square;
set(gca,'fontsize',Font_size);
xlabel('t');
title('Energy');

subplot(2,3,4);
mesh(x,y,wrapToPi(angle(psiL_temp)));
axis square;
set(gca,'fontsize',Font_size,'zlim',[-pi, pi]);
view(0,90);
colorbar;
zoom(zoom_value);
xlabel('x');
ylabel('y');
title('\phi_L');

subplot(2,3,5);
mesh(x,y,wrapToPi(angle(psiR_temp)));
axis square;
set(gca,'fontsize',Font_size,'zlim',[-pi, pi]);
view(0,90);
colorbar;
zoom(zoom_value);
xlabel('x');
ylabel('y');
title('\phi_R');

%calculate velocity field [vL_x, vL_y], [vR_x, vR_y]
hspace=Ori_Pump_Parameter.hspace;
mod_psiL=abs(psiL_temp);
mod_psiR=abs(psiR_temp);

[psix, psiy]=gradient(psiL_temp,hspace);
[psix_c, psiy_c]=gradient(conj(psiL_temp),hspace);

vL_x=(conj(psiL_temp).*psix-psiL_temp.*psix_c)./(1i*mod_psiL.^2);
vL_y=(conj(psiL_temp).*psiy-psiL_temp.*psiy_c)./(1i*mod_psiL.^2);

[psix, psiy]=gradient(psiR_temp,hspace);
[psix_c, psiy_c]=gradient(conj(psiR_temp),hspace);

vR_x=(conj(psiR_temp).*psix-psiR_temp.*psix_c)./(1i*mod_psiR.^2);
vR_y=(conj(psiR_temp).*psiy-psiR_temp.*psiy_c)./(1i*mod_psiR.^2);

%calculate vortexity

[curlLz,cavL]= curl(x,y,vL_x,vL_y);
[curlRz,cavR]= curl(x,y,vR_x,vR_y); %curlRz vorticity z component, cavR angular velocity

%plot vortexity of R component

x_b=10;%plot range for the untrimmed
N1=floor((XYmax-x_b)/hspace);
N2=floor((XYmax+x_b)/hspace);

subplot(2,3,6);
mesh(x(N1:N2,N1:N2),y(N1:N2,N1:N2),curlRz(N1:N2,N1:N2));
view(0,30);
axis square;
colorbar;
set(gca,'fontsize',Font_size,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
xlabel('x');
ylabel('y')
title('Vortexity of \psi_R');

print(gcf,'-dpng','-zbuffer','-r200','phase_and_vorticity.png');









