fig_E=figure('position',[pos_1(1)+Windows_shift*6   333   windows_width+200   windows_hight+100],'renderer','painter','paperpositionmode','auto');
[AX,H1,H2]=plotyy(Time,Energy,Time,real(Lz));
set(AX(1),'fontsize',Font_size,'fontname',Font_name);
set(AX(1),'xlim',[0,max(Time)]);
set(AX(1),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(1),'ylim',[1 7],'ytick',1:7);
set(H1,'linewidth',3);
xlabel('t');
ylabel(AX(1),'Energy');
box(AX(1),'off');

set(AX(2),'fontsize',Font_size,'fontname',Font_name);
set(AX(2),'xlim',[0,max(Time)]);
set(AX(2),'xtick',[]);
set(AX(2),'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
set(AX(2),'ylim',[10 80]);
set(AX(2),'ytick',10:20:80);
set(H2,'linewidth',3,'linestyle','--');
ylabel(AX(2),'L_z');
box(AX(2),'off');

linkaxes(AX,'x');
set(AX(2), 'XTickLabel','','XAxisLocation','Top') 

ann_position=[ 0.7531    0.8267    0.0800    0.0800];
h_a=annotation('textbox', ann_position,'String','(b)');
set(h_a,'BackgroundColor',[1 1 1],'Color','k','FontName',Font_name,'fontsize',Font_size+4,...
    'HorizontalAlignment','center','VerticalAlignment','top','LineStyle','none','FitHeightToText','on');
set(h_a,'position',ann_position);