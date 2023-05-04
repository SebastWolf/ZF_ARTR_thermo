%% =========================== Supplementary Figure 7A-B ============================================
%Kullback-Leibler divergence between the experimental and the Langevin distributions as a function of A for two data sets.

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];
load([pwd filesep 'Data availability' filesep 'Mean_field_Parameters.mat'])

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

id_fish = 10

dpath = hpath_star(fish_num(id_fish));
clear files_T
files = dir(dpath);

for t_files = [3 5];
    for idT = 1:numel(T)
        if isempty(strfind(files(t_files).name,num2str(T(idT))))==0;
            files_T(t_files) = T(idT);
        end
    end
    idT = find(files_T(t_files)==T);
    
    load([files(t_files).folder filesep files(t_files).name]);
    
    mleft_data = Dinference_corr.mleft_data;
    mright_data = Dinference_corr.mright_data;
    
    
    L_reg = Dinference_corr.L_reg;
    R_reg = Dinference_corr.R_reg;
    N_L = length(L_reg);
    N_R = length(R_reg);
    
    
    h_L = mean(Dinference_corr.h(1:length(L_reg)));
    h_R = mean(Dinference_corr.h(length(L_reg)+1:length(L_reg)+length(R_reg)));
    
    std_h_L = std(Dinference_corr.h(1:length(L_reg)));
    std_h_R = std(Dinference_corr.h(length(L_reg)+1:length(L_reg)+length(R_reg)));
    
    J_inter = Dinference_corr.J(1:N_L,N_L+1:N_L+N_R);
    J_left = Dinference_corr.J(1:N_L,1:N_L);
    J_right = Dinference_corr.J(N_L+1:N_L+N_R,N_L+1:N_L+N_R) ;
    
    for i = 1:size(J_left,1); J_left(i,i) = NaN;end
    for i = 1:size(J_right,1); J_right(i,i) = NaN;end
    
    J_L = nanmean(J_left(:))*(length(L_reg));
    J_R = nanmean(J_right(:))*(length(R_reg));
    I_I = nanmean(J_inter(:))*sqrt((length(R_reg))*(length(L_reg)));
    U = triu(J_left)
    U(U==0) = NaN;
    
    %     % compute landscape
    sig_h_R = std_h_R;
    sig_h_L = std_h_L;
    
    M_h_R = h_R
    M_J_R = nanmean(J_right(:))
    
    M_h_L = h_L;
    M_J_L = nanmean(J_left(:))
    
    M_I_I = nanmean(J_inter(:));
    
    J_J_L = nanmean(J_left(:))*(length(L_reg));
    J_J_R = nanmean(J_right(:))*(length(R_reg));
    I_I = nanmean(J_inter(:))*sqrt((length(R_reg))*(length(L_reg)));
    
    % activity
    N_L_keep = N_L
    N_R_keep = N_R
    
    clear cc_list K_list factor_list
    
    factor_list = [2:0.1:max(N_L_keep,N_R_keep)/4];
    
    for alpha = 1:length(factor_list);
        
        factor = factor_list(alpha);
        %factor = 30
        N_L_noise = N_L_keep/factor;
        N_R_noise = N_R_keep/factor;
        
        % generator of activity
        
        %         K_L = 20;
        %         K_R = 20;
        tau = 1;
        dt = 0.001;
        
        thr_lim = 1e-4;
        
        sigma_L = sqrt(tau/dt)*sqrt(2);
        sigma_R = sqrt(tau/dt)*sqrt(2);
        T_time = 200000;
        mr = thr_lim*ones(1,T_time);
        ml = thr_lim*ones(1,T_time);
        
        N_L = N_L_noise;
        N_R = N_R_noise;
        
        
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
        
        %data
        nbins = 10;
        pseudocount = 1;
        ebins = linspace(0, 1, nbins + 1);    % bin edges
        
        % Get counts in each bins
        countdata = histcounts2(mleft_data, mright_data, ebins, ebins);
        countisin = histcounts2(ml, mr, ebins, ebins);
        
        % Add pseudocounts
        countdata = countdata + pseudocount;
        countisin = countisin + pseudocount;
        
        % Normalize to get proba
        pdfdata = countdata./sum(countdata, 'all'); % proba
        pdfisin = countisin./sum(countisin, 'all');
        
        % Linearize
        pdfdata = pdfdata(:);
        pdfisin = pdfisin(:);
        
        % KL divergence
        K_list(alpha) = sum( pdfdata.*log10(pdfdata./pdfisin) );
        cc = corrcoef(countisin(:),countdata(:));
        cc_list(alpha) = cc(1,2);
        
    end
    
    
    figure;
    plot(factor_list,smooth(K_list,10),'k','LineWidth',3)
    set(gca,'FontSize',18)
    xlabel('A')
    ylabel('DKL')
    ylim([0 0.6])
    box off
    
    
end
