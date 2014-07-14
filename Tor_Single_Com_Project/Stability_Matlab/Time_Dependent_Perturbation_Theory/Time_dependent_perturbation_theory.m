
%Time dependent perturbation theory

Windows_width=800;
Windows_hight=Windows_width;

Font_size=20;

print_figure=1;

r0=16*sqrt(5/2);

m=60;
k=65;

o1=1.9445;
o2=0.1233;


t_series=0.1:0.01:30;

fig_Et=figure('position',[500 300 Windows_width Windows_hight],'renderer','painters','paperpositionmode','auto');

P=(1+exp(2*t_series*o2)-2*exp(t_series*o2).*cos(t_series*((m^2-(m-k)^2)/(2*r0^2)+o1)))./t_series;

plot(t_series,P,'LineWidth',1);
axis square
set(gca,'fontsize',Font_size);
xlabel('t');
ylabel('bare transition rate P/t');
title(sprintf('r0=%.0f, m0=%d, k=%.0f, omega1=%.2f, omega2=%.2f',r0,m,k,o1,o2));


if print_figure~=0
    
   print(fig_Et,'-painters','-dpng','-r250',sprintf('Transition_rate_r0=%.0f, m0=%d, k=%.0f, omega1=%.2f, omega2=%.2f.png',r0,m,k,o1,o2))
   
   close(fig_Et);
    
end

















