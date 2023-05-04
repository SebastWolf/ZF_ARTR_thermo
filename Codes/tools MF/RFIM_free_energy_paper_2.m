function F =  RFIM_free_energy_paper_2(u,v,N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I);
%M_I_I = moyenne sur popu croisées
%M_h_L = moyenne sur popu left
%M_h_R = moyenne sur popu right


H_L = M_h_L-M_J_L/2;
H_R = M_h_R-M_J_R/2;

F = -N_L*J_J_L/2*u^2-N_R*J_J_R/2*v^2-I_I*sqrt(N_L*N_R)*u*v-N_L*H_L*u-N_R*H_R*v+N_L*(u*log(u)+(1-u)*log(1-u))+N_R*(v*log(v)+(1-v)*log(1-v));


end