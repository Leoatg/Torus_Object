

fig_E=figure('position',[50 80 1800 900],'PaperPositionMode','auto','renderer','opengl','Visible','on');


subplot(2,2,1);
mesh(x,y,abs(psi).^2)
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
set(gca,'fontname',Font_name);
xlabel('x');
ylabel('y');
view(0,90);
zoom(Psi.zoom_factor);
colormap jet;
colorbar;
title('|\psi|^2');

subplot(2,2,3)
mesh(x,y,angle(psi));
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax],'zLim',[-pi pi]);
set(gca,'fontname',Font_name);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('\phi');

subplot(2,2,2)
plot(Time,Energy_unmal,'linewidth',2);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
pos_tem2=pos_tem;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size);
set(gca,'fontname',Font_name);
xlabel('t');
%ylabel('y');
title('Energy (Unnormalized)');

subplot(2,2,4);
mesh(x,y,n);
pos_tem=get(gca,'position');
pos_tem(1)=0.3692;
set(gca,'position',pos_tem);
axis square;
set(gca,'fontsize',Psi.font_size,'xLim',[-XYmax XYmax],'yLim',[-XYmax XYmax]);
set(gca,'fontname',Font_name);
view(0,90);
zoom(Psi.zoom_factor);
colorbar;
xlabel('x');
ylabel('y');
title('n');

%right bottom corner
axes('position',[ 0.6420    0.1256    0.3347    0.3412]);
plot(Time,Energy,'linewidth',2);
pos_tem(1)=0.5859;
set(gca,'position',pos_tem);
set(gca,'fontsize',Psi.font_size);
set(gca,'fontname',Font_name);
if max(num2str(get(gca,'YTick')))<1
    set(gca,'YTickLabel',num2str(get(gca,'YTick')','%.2f'));
end
axis square;
xlabel('t');
title('Energy (normalized)');

%right upper corner
pos_tem2(1)=0.5859;
axes('position',pos_tem2);
plot(Time,real(Lz),'linewidth',2);
set(gca,'fontsize',Psi.font_size);
set(gca,'fontname',Font_name);
axis square;
xlabel('t');
title(sprintf('L_z (normalized)  m_0=%d;  L_z final=%.1f',Psi.Ini_m,real(Lz(length(Lz)))));

