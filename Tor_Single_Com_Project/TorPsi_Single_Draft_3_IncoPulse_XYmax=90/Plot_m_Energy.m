
%Plot m vs Energy

m_series=0:50;
N_m=length(m_series);

Energy_m_nor=zeros(1,N_m);
Energy_m_unnor=Energy_m_nor;
Psi_int_m=Energy_m_nor;

for it=1:N_m
    
    m=m_series(it);
    
    file_name=sprintf('FinalState_TorPsi_IncoPulse_GmaRr=1.5_m0=%.1f_Pmp0=2.50_l=5_p=0_w0=28_XYmax=90_N=1024_T=150.mat',m);
    
    load(file_name,'Time','Energy','Psi_int');
    
    N_t=length(Time);
    
    Energy_m_nor(it)=Energy(N_t);
    
    Energy_m_unnor(it)=Energy(N_t)*Psi_int(N_t);
    
    Psi_int_m(it)=Psi_int(N_t);
end

figure;
plot(m_series,Energy_m_nor,'linewidth',2);
xlabel('m');
title('Energy (normalized)');

figure;
plot(m_series,Energy_m_unnor,'linewidth',2);
xlabel('m');
title('Energy (Unnormalized)');

figure;
plot(m_series,Psi_int_m,'linewidth',2);
xlabel('m');
title('Psi int');









