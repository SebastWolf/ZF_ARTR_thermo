function DM_U_free_energy =  RFIM_free_energy_DM_U_paper_2(u,v,N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I);

%I_I = M_I_I*sqrt(N_L*N_R);
%J_J_L = N_L*M_J_L;
%J_J_R = N_R*M_J_R;

H_L = M_h_L-M_J_L/2;
H_R = M_h_R-M_J_R/2;
DM_U_free_energy = (-N_L*(J_J_L*u+H_L)+I_I*sqrt(N_L*N_R)*v+N_L*log(u/(1-u)));
%DM_U_free_energy = (-(J_J_L*u+H_L)+I_I*sqrt(N_L*N_R)/N_L*v+log(u/(1-u)));

end