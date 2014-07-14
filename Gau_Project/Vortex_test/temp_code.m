

%temp



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









