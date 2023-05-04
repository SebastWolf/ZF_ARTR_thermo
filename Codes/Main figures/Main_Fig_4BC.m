%% ===========================  Figure 4B-C ============================================
% Examples of simulated mL and mR signals of the mean-field
% dynamical equations for three sets of parameters corresponding to three bath temperatures. 
% Free-energy landscapes computed with the mean-field model, in the (mL,mR) plane for three different 
% parameter sets.

clear all
clc
close all
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
load([pwd '/Data availability/Mean_field_Parameters.mat'])
hpath_star = @(fish_num) [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J/T*_Fish ' num2str(fish_num) '.mat'];
fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

id_fish = 4

dpath = hpath_star(fish_num(id_fish));
clear files_T
files = dir(dpath);

for t_files = 1:length(files)
    for idT = 1:numel(T)
        if isempty(strfind(files(t_files).name,num2str(T(idT))))==0;
            files_T(t_files) = T(idT);
        end
    end
    idT = find(files_T(t_files)==T);
    
    line_selected = find(table_mean_field(:,1)==T(idT) & table_mean_field(:,2)==fish_num(id_fish));
    
    load([files(t_files).folder filesep files(t_files).name]);
    
    L_reg = Dinference_corr.L_reg;
    R_reg = Dinference_corr.R_reg;
    N_L = length(L_reg);
    N_R = length(R_reg);
    J_inter = Dinference_corr.J(1:N_L,N_L+1:N_L+N_R);
    J_left = Dinference_corr.J(1:N_L,1:N_L);
    J_right = Dinference_corr.J(N_L+1:N_L+N_R,N_L+1:N_L+N_R) ;
    
    %     % compute landscape
    
    M_h_R = table_mean_field(line_selected,7);
    M_h_L = table_mean_field(line_selected,6);;
    
    M_J_R = nanmean(J_right(:))
    M_J_L = nanmean(J_left(:))
    M_I_I = nanmean(J_inter(:));
    
    J_J_L = table_mean_field(line_selected,3);
    J_J_R = table_mean_field(line_selected,4);
    I_I = table_mean_field(line_selected,5);
    
    % activity
        
    N_L = table_mean_field(line_selected,8);
    N_R = table_mean_field(line_selected,9);
    
    u_list = [1e-4:0.005:1-1e-4];
    v_list = [1e-4:0.005:1-1e-4];
    
    
    clear free_energy_all
    
    for i = 1:length(u_list);
        
        for j = 1:length(v_list);
            u = u_list(i);
            v = v_list(j);
            free_energy_all(i,j) = RFIM_free_energy_paper_2(u,v,N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I);
        end
    end
    
    figure
    sss = surf(free_energy_all)
    sss.EdgeColor = 'none';
    title(['T = ' num2str(files_T(t_files))])
    set(gca,'FontSize',18)
    
    
    tau = 1;
    dt = 0.001;
    
    thr_lim = 1e-4;
    
    sigma_L = sqrt(tau/dt)*sqrt(2)
    sigma_R = sqrt(tau/dt)*sqrt(2);
    T_time = 10000;
    mr = thr_lim*ones(1,T_time);
    ml = thr_lim*ones(1,T_time);
    
    
    for i = 1:T_time-1
        
        %ml
        xi_L = sigma_L*randn(1);
        xi_L_note(i) = xi_L;
        ml(i+1) = ml(i)-dt*(1/tau)*(RFIM_free_energy_DM_U_paper_2(ml(i),mr(i),N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I)+xi_L);
        if ml(i+1)<=thr_lim; ml(i+1) = thr_lim; end
        if ml(i+1)>=1-thr_lim; ml(i+1) = 1-thr_lim; end
        
        %mr
        xi_R = sigma_R*randn(1);
        mr(i+1) = mr(i)-dt*(1/tau)*(RFIM_free_energy_DM_V_paper_2(ml(i),mr(i),N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I)+xi_R);
        if mr(i+1)<=thr_lim; mr(i+1) = thr_lim; end
        if mr(i+1)>=1-thr_lim; mr(i+1) = 1-thr_lim; end
        
    end
    
    
    vect_time = dt*[1:T_time];
    figure('Name',files(t_files).name)
    plot(vect_time,smooth(mr,5)); hold on; plot(vect_time,smooth(ml,5))
    title(['T = ' num2str(files_T(t_files))])
    set(gca,'FontSize',14)
    
    
end
