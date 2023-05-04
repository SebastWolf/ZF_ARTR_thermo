%% ================================ Supplementary Figure 85BC ============================================

clear all
clc
close all


pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];

T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

fish_out = [];
temp_out = [];

h_L_ALL = [];
h_R_ALL = [];
J_LR_ALL = [];

pcount = 0;

for id_fish = 1:numel(fish_num)
    
    dpath = hpath_star(fish_num(id_fish));
    files = dir(dpath);
    
    clear files_T
    for t_files = 1:length(files)
        for idT = 1:numel(T)
            if isempty(strfind(files(t_files).name,num2str(T(idT))))==0;
                files_T(t_files) = T(idT);
            end
        end
        idT = find(files_T(t_files)==T);
        
        id_fish;
        t_files
        
            
            load([files(t_files).folder filesep files(t_files).name]);
            
            h_ = Dinference_corr.h;
            h_L = h_(Dinference_corr.L_reg);
            h_R = h_(Dinference_corr.R_reg);
            
            J_ = Dinference_corr.J;
            for i=1:size(J_,1); J_(i,i) = NaN; end
            
            J_L = J_(Dinference_corr.L_reg,Dinference_corr.L_reg);
            J_R = J_(Dinference_corr.R_reg,Dinference_corr.R_reg);
            J_LR = J_(Dinference_corr.L_reg,Dinference_corr.R_reg);
            
            h_L_ALL = [h_L_ALL h_L(:)'];
            h_R_ALL = [h_R_ALL h_R(:)'];
     
            pcount = pcount+1;
            m_h(pcount,1) = nanmean(h_L);
            m_h(pcount,2) = nanmean(h_R);
            
            m_std(pcount,1) = nanstd(h_L);
            m_std(pcount,2) = nanstd(h_R);
            
    end
    
end


%%

figure; boxplot(m_h)
set(gca,'FontSize',18)
box off

figure; boxplot(m_std)
set(gca,'FontSize',18)
box off

%%

mean([h_L_ALL h_R_ALL])
std([h_L_ALL h_R_ALL])