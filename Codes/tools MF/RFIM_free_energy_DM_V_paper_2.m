function DM_V_free_energy =  RFIM_free_energy_DM_V_paper_2(u,v,N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I);

%I_I = M_I_I*sqrt(N_L*N_R);
%J_J_L = N_L*M_J_L;
%J_J_R = N_R*M_J_R;

H_L = M_h_L-M_J_L/2;
H_R = M_h_R-M_J_R/2;

DM_V_free_energy = (-N_R*(J_J_R*v+H_R)+I_I*sqrt(N_L*N_R)*u+N_R*log(v/(1-v)));
%DM_V_free_energy = (-(J_J_R*v+H_R)+I_I*sqrt(N_L*N_R)/N_R*u+log(v/(1-v)));

end