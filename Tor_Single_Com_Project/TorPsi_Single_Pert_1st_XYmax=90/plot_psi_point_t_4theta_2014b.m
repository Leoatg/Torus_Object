
%plot psi point

C=Control_Parmeters;

windows_width=500;
windows_hight=500;

Font_size=20;

N_pp=length(C.psi_points_theta);

print_fig=0;

fig_psi_point=figure('position',[200 107 1017 985],'renderer','painters','paperpositionmode','auto'); 


for it=1:N_pp
    
   subplot(2,2,it);
   plot(Time,abs(record_psi_points(it,:)).^2,'LineWidth',1.5);
   axis square;
   set(gca,'fontsize',Font_size);
   set(gca,'xlim',[0 max(Time)]);
   xlabel('t');
   %set(gca,'ytick',floor(get(gca,'ytick')));
   title(['|\psi|^2 at \theta=',num2str(C.psi_points_theta(it),'%.2f')]);
    
end

if print_fig~=0
    
   print(fig_psi_point,'-painters','-dpng','-r200',sprintf('m0=%d,psi_4points.png',round(real(Lz(1)))));
    close(fig_psi_point);
end








