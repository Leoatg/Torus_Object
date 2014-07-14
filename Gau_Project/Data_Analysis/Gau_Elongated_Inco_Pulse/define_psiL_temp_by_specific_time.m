

%define psiL_temp/psiR_temp
try 
prompt='\n input psiL_temp/psiR_temp t:  ';
psi_temp_t = input(prompt);
end

addpath('D:\Torus_Pulse_Object');

if ~exist('record_psi_t','var')
    fprintf('Error: No psi recorded\n');
    return;
elseif psi_temp_t<min(record_psi_t) || psi_temp_t>max(record_psi_t)
    fprintf('Error: psiL_temp/psiR_temp time error\n')
    fprintf('time value should between:%.2f - %.2f\n',min(record_psi_t),max(record_psi_t));
    return;
end

psi_temp_num=length(min(record_psi_t):(record_psi_t(2)-record_psi_t(1)):psi_temp_t);

psiL_temp=psiL_trim(:,:,psi_temp_num);psiR_temp=psiR_trim(:,:,psi_temp_num);

N=Ori_Pump_Parameter.N;
XYmax=Ori_Pump_Parameter.XYmax;
hspace=Ori_Pump_Parameter.hspace;
dt=Ori_Pump_Parameter.dt;

[x,y]=meshgrid(-Ori_Pump_Parameter.XYmax:Ori_Pump_Parameter.hspace:Ori_Pump_Parameter.XYmax);

x_axis=-Ori_Pump_Parameter.XYmax:Ori_Pump_Parameter.hspace:Ori_Pump_Parameter.XYmax;


fprintf('psiL_temp=psiL_trim(:,:,%d);psiR_temp=psiR_trim(:,:,%d);\n',psi_temp_num,psi_temp_num)

fprintf('\n');

clear psi_temp_t psi_temp_num prompt N_psiL_temp