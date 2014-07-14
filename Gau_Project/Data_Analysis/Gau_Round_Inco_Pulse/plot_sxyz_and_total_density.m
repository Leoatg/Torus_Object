
%plot sxyz int and Total density

Font_size=16;

figure('position',[100 100 1200 800],'renderer','painter','paperpositionmode','auto');

subplot(2,2,1);
[AX,H1,H2]=plotyy(Time,sx_int,Time,mod_psiLR_square);
set(AX(1),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1],'xtick',[]);
set(H1,'linewidth',1.5);
set(H2,'linewidth',1.5);
xlabel(AX(1),'t');
ylabel(AX(1),'sx int');
ylabel(AX(2),'Total Density');

subplot(2,2,2);
[AX,H1,H2]=plotyy(Time,sy_int,Time,mod_psiLR_square);
set(AX(1),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1],'xtick',[]);
set(H1,'linewidth',1.5);
set(H2,'linewidth',1.5);
xlabel(AX(1),'t');
ylabel(AX(1),'sy int');
ylabel(AX(2),'Total Density');

subplot(2,2,3);
[AX,H1,H2]=plotyy(Time,sz_int,Time,mod_psiLR_square);
set(AX(1),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1],'xtick',[]);
set(H1,'linewidth',1.5);
set(H2,'linewidth',1.5);
xlabel(AX(1),'t');
ylabel(AX(1),'sz int');
ylabel(AX(2),'Total Density');

subplot(2,2,4);
[AX,H1,H2]=plotyy(Time,Energy,Time,mod_psiLR_square);
set(AX(1),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'fontsize',Font_size,'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1],'xtick',[]);
set(H1,'linewidth',1.5);
set(H2,'linewidth',1.5);
xlabel(AX(1),'t');
ylabel(AX(1),'Energy');
ylabel(AX(2),'Total Density');

print(gcf,'-painter','-dpng','-r200','sxyz_and_TotalDensity.png');







