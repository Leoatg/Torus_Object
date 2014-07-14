

%Build mesh

N=Ori_Pump_Parameter.N;
XYmax=Ori_Pump_Parameter.XYmax;
hspace=2*XYmax/(N-1);

[x,y]=meshgrid(-XYmax:hspace:XYmax,-XYmax:hspace:XYmax);

psi=Finial_State(:,:,1);
psi_mode=abs(psi);
psi_angle=angle(psi);
theta=atan2(y,x);