
%Pattern analysis

clear Pat I_std phi_arc phi_mod

Pat=Pat_test;

Pat=abs(Pat).^2;

Pat=Pat/max(max(Pat)); %intensity normalize to 0-1

Np=length(Pat);

N_wx=14;

N_wy=N_wx;

for it=N_wx+1:(Np-N_wx-1)
    for jt=N_wy+1:(Np-N_wy-1)
        
        Nbx=(jt-N_wx);
        Nex=(jt+N_wx);
        Nby=(it-N_wx);
        Ney=(it+N_wx);
        
        Max=max(max(Pat(Nby:Ney,Nbx:Nex)));
        Min=min(min(Pat(Nby:Ney,Nbx:Nex)));
    
        I_std(it-N_wy,jt-N_wx)=(Pat(it,jt)-Min)/(Max-Min);
    
    end
end

x_t=x(N_wx+1:(Np-N_wx-1),N_wx+1:(Np-N_wx-1));
y_t=y(N_wx+1:(Np-N_wx-1),N_wx+1:(Np-N_wx-1));

phi_arc=acos(2*I_std-1);

%turn phi_arc to prime value

N_arc=length(phi_arc);
phi_mod=phi_arc;

delta=2*pi/N_arc;
%delta=0.1*pi;
 delta_up=pi-0.2;
 delta_low=0.9;%delta;

upper_end_flag=0;
lower_end_flag=0;
phase_inversion_flag=0;

for it=1:(N_arc) %for each line
%     if (it<100) || (it> (N_arc-100))
%     delta_up=min(findpeaks(phi_arc(it,:)))*0.9;
%     delta_low=min(phi_arc(it,:))+delta;
%     else
%         delta_up=pi*0.8;
%         delta_low=0.6;
%     end
    
    upper_end_flag=0;
    lower_end_flag=0;
    phase_inversion_flag=0;
    
    for jt=1:(N_arc-1)
        
        if phi_arc(it,jt)<delta_low
            upper_end_flag=0;
            lower_end_flag=1;
        elseif phi_arc(it,jt)>delta_up
            upper_end_flag=1;
            lower_end_flag=0;
        end
        
        if upper_end_flag==1
            if (phi_arc(it,jt+1)-phi_arc(it,jt))<0
                phase_inversion_flag=1;
            end
        end
        if lower_end_flag==1
            if (phi_arc(it,jt+1)-phi_arc(it,jt))>0
                phase_inversion_flag=0;
            end
            
        end
        
        if phase_inversion_flag==0
            phi_mod(it,jt)=phi_arc(it,jt);
        else
            phi_mod(it,jt)=2*pi-phi_arc(it,jt);
        end
        
    end
    
end

fig_Pat_ana=figure('renderer','zbuffer','position',[100+Box_width 504 Box_width+50 Box_width ],'paperpositionmode','auto');
mesh(x_t,y_t,phi_mod);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
xlabel('x');
ylabel('y');
view(0,90);
colormap gray;
colorbar;
title('Pattern phase prime value');

% phi_mod_un=unwrap(unwrap(wrapToPi(phi_mod),[],1),[],2);
% 
% fig_phi_mod_un=figure('renderer','zbuffer','position',[100+Box_width 504 Box_width+50 Box_width ],'paperpositionmode','auto');
% mesh(x_t,y_t,phi_mod_un);
% axis square;
% set(gca,'fontsize',Font_size,'xlim',[-XYmax XYmax],'ylim',[-XYmax XYmax]);
% xlabel('x');
% ylabel('y');
% view(0,90);
% colormap gray;
% colorbar;
% title('Pattern phase unwrapped');


%trim phase prime value

x_b=2.5;



N1=floor((XYmax-x_b)/hspace);
N2=floor((XYmax+x_b)/hspace);

phi_mod_t=phi_mod(N1:N2,N1:N2);

x_t_2=x_t(N1:N2,N1:N2);
y_t_2=y_t(N1:N2,N1:N2);
XYmax_b=max(max(x_t_2));

phi_mod_t=unwrap(unwrap(phi_mod_t,[],1),[],2);

fig_phi_mod_trim_un=figure('renderer','zbuffer','position',[100+Box_width 504 Box_width+50 Box_width ],'paperpositionmode','auto');
mesh(x_t_2,y_t_2,phi_mod_t);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax_b XYmax_b],'ylim',[-XYmax_b XYmax_b]);
xlabel('x');
ylabel('y');
view(0,90);
colormap gray;
colorbar;
title('Pattern phase unwrapped');


%subtract background phase

% fit at y=-1
Ny=floor((XYmax_b-1)/hspace);
x_axis=x_t_2(Ny,:);
Tar_line=phi_mod_t(Ny,:);
f=fit(x_axis',Tar_line','poly1');
%Coe=coeffnames(f);
Coe_v=coeffvalues(f);
BackGround=Coe_v(1)*x_t_2+Coe_v(2);

fig_background=figure('renderer','zbuffer','position',[200 161 558 480],'paperpositionmode','auto');
mesh(x_t_2,y_t_2,BackGround);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax_b XYmax_b],'ylim',[-XYmax_b XYmax_b]);
xlabel('x');
ylabel('y');
view(0,90);
colormap gray;
colorbar;
title('Background Phase')

Vortex_phase=phi_mod_t-BackGround;

fig_vortexphase=figure('renderer','zbuffer','position',[300 361 558 480],'paperpositionmode','auto');
mesh(x_t_2,y_t_2,Vortex_phase);
axis square;
set(gca,'fontsize',Font_size,'xlim',[-XYmax_b XYmax_b],'ylim',[-XYmax_b XYmax_b]);
xlabel('x');
ylabel('y');
view(0,90);
colormap gray;
colorbar;
title('Vortex Phase');





% print(fig_Pat_ana,'-zbuffer','-dpng','Pattern phase prime value, view from top.png')
% 
% print(fig_phi_mod_trim_un,'-zbuffer','-dpng','Pattern phase unwrapped, view from top.png')


















