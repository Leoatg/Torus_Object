

%plot phase

Font_size=16;
%plot phase
figure('position',[100 100 1500 800],'renderer','zbuffer','paperpositionmode','auto');

phiR=wrapToPi(angle(psiR_temp));
subplot(1,3,1);
mesh(x,y,phiR);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-x_b x_b],'ylim',[-x_b,x_b],'zlim',[-pi pi]);
xlabel('x');
ylabel('y');
view(0,90);
colorbar;
title('\phi_R wrapped');

unwrapped_phase=unwrap(unwrap(phiR,[],1),[],2);
unwrapped_min=min(min(unwrapped_phase));
subplot(1,3,2);
mesh(x,y,unwrapped_phase-phiR_Steady);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-x_b x_b],'ylim',[-x_b,x_b]);
%colorbar;
%caxis([0 3*pi]);
view(60,8);
xlabel('x');
ylabel('y');
title('unwrapped \phi_R - unwrapped Steady State phase');

subplot(1,3,3);
mesh(x,y,unwrapped_phase-unwrapped_min);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-x_b x_b],'ylim',[-x_b,x_b]);
xlabel('x');
ylabel('y');
%caxis([0 3*pi]);
colorbar;
view(0,90);
title('unwrapped \phi_R, viewed from top');

print(gcf,'-zbuffer','-dpng','-r150','phase.png');






