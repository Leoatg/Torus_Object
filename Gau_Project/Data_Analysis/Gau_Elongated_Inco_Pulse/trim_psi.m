

%trim psi as indicated x_b

x_b=16;

N1=floor((XYmax-x_b)/hspace);
N2=floor((XYmax+x_b)/hspace);

psiL_t=psiL_temp(N1:N2,N1:N2);
psiR_t=psiR_temp(N1:N2,N1:N2);
x_t=x(N1:N2,N1:N2);
y_t=y(N1:N2,N1:N2);

clear psiL_temp psiR_temp x y
psiL_temp=psiL_t;
psiR_temp=psiR_t;
x=x_t;
y=y_t;



%unwrap the phase

% N_b=N2-N1+1;
% 
% phase1=zeros(1,N_b^2);
% phase_u=phase1; %unwraped 1D phase
% D=zeros(1,length(phase1)-1);
% Delta=D;
% phase_unwrap=zeros(N_b,N_b);
% c=1;
% for i=1:N_b %row sweep
%     
%     if mod(i,2)~=0
%     
%         for j=1:N_b
%             phase1(c)=angle(psiR_t(i,j));
%             c=c+1;
%         end
%         
%     else 
%         for j=N_b:-1:1
%             phase1(c)=angle(psiR_t(i,j));
%             c=c+1;
%         end
%     end
%     
% end
% 
% for i=1:length(D)
%     D(i)=phase1(i+1)-phase1(i);
% end
% 
% for i=1:length(D)
%     Delta(i)=atan(sin(D(i))/cos(D(i)));
% end
% 
% for i=1:(length(D)-1)
%     phase_u(i+1)=phase_u(i)+Delta(i+1);
% end
% 
% c=1;
% for i=1:N_b %row sweep
%     
%     if mod(i,2)~=0
%     
%         for j=1:N_b
%             phase_unwrap(i,j)=phase_u(c);
%             c=c+1;
%         end
%         
%     else 
%         for j=N_b:-1:1
%             phase_unwrap(i,j)=phase_u(c);
%             c=c+1;
%         end
%     end
% end
% 
% 
% f_psi1=figure('renderer','zbuffer','paperpositionmode','auto');
% mesh(x_t,y_t,abs(psiR_t).^2);
% set(gca,'fontsize',16,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
% view(0,90);
% colorbar;
% xlabel('x');
% ylabel('y');
% title('|\psi_r|^2');
% 
% f_psi2=figure('renderer','zbuffer','paperpositionmode','auto');
% mesh(x_t,y_t,abs(psiR_t).^2);
% set(gca,'fontsize',16,'xlim',[-x_b,x_b],'ylim',[-x_b,x_b]);
% %view(0,90);
% colorbar;
% xlabel('x');
% ylabel('y');
% title('|\psi_+|^2');
% 
% x_axis=-x_b:hspace:x_b;
% f_psi3=figure('paperpositionmode','auto');
% [AX,H1,H2]=plotyy(x_axis,abs(psiR_t(floor(N_b/2)+1,:)).^2,x_axis,phase_unwrap(floor(N_b/2)+1,:));
% set(AX(1),'fontsize',16,'xlim',[-x_b,x_b],'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
% set(AX(2),'fontsize',16,'xlim',[-x_b,x_b],'xtick',[],'DataAspectRatioMode','auto','PlotBoxAspectRatio',[1 1 1]);
% set(H1,'linewidth',1.5);
% set(H2,'linewidth',1.5);
% h=legend('|\psi_+|^2','\phi_+');
% set(h,'Orientation','horizontal','box','off','position',[0.3351    0.92    0.3155    0.0858]);
% 
% print(f_psi1,'-zbuffer','-r130','-depsc','fpsi1.eps');
% print(f_psi2,'-zbuffer','-r130','-depsc','fpsi2.eps');
% print(f_psi3,'-painter','-r200','-depsc','fpsi3.eps');









