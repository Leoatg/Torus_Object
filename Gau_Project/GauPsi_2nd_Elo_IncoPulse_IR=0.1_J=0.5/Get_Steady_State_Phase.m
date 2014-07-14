
%Get Steady State phase

addpath('D:\Torus_Pulse_Object');

load('Gau_Standard_Steady_State_2nd_Elo_m0_0.mat');

psiL_temp=Finial_State(:,:,1);
psiR_temp=Finial_State(:,:,2);

N=Ori_Pump_Parameter.N;
XYmax=Ori_Pump_Parameter.XYmax;
hspace=Ori_Pump_Parameter.hspace;
dt=Ori_Pump_Parameter.dt;

[x,y]=meshgrid(-Ori_Pump_Parameter.XYmax:Ori_Pump_Parameter.hspace:Ori_Pump_Parameter.XYmax);

x_axis=-Ori_Pump_Parameter.XYmax:Ori_Pump_Parameter.hspace:Ori_Pump_Parameter.XYmax;


%trim

x_b=10;

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

phiR=wrapToPi(angle(psiR_temp));
phiR_Steady=unwrap(unwrap(phiR,[],1),[],2);

fprintf('\n phiR_Steady\n');


