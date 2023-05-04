%% ================================ SI Figure 6 ============================================

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];
pcount_total = 0;

for id_fish = 1:numel(fish_num)
    clear J_ h_ J_sum_
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
        
        load([files(t_files).folder filesep files(t_files).name]);
        
        J_(:,:,t_files) = Dinference_corr.J ;
        h_(:,t_files) = Dinference_corr.h ;
        
        J_sum_(:,t_files) = sum(Dinference_corr.J,2)';
        
    end
    
    for t_files1 = 1:length(files)
        for t_files2 = t_files1+1:length(files)
            
            X_J = squeeze(J_(:,:,t_files1));
            Y_J = squeeze(J_(:,:,t_files2));
            X_h = squeeze(h_(:,t_files1));
            Y_h = squeeze(h_(:,t_files2)); 
            X_J_sum = squeeze(J_sum_(:,t_files1))+X_h;
            Y_J_sum = squeeze(J_sum_(:,t_files2))+Y_h;
           
            pcount_total = pcount_total + 1
                                                            
            R = corrcoef([X_J(:),Y_J(:)]);
            R_J_keep(pcount_total) = R(1,2);
            
            R = corrcoef([X_J_sum(:),Y_J_sum(:)]);
            R_J_sum_keep(pcount_total) = R(1,2);
            
            R = corrcoef([X_h(:),Y_h(:)]);
            R_h_keep(pcount_total) = R(1,2);
            
        end
    end
end

%%

figure;
boxplot([R_h_keep.^2' R_J_keep.^2'])
%ylabel('Training vs test (%)')
set(gca,'FontSize',18)
box off
xticks([1 2])
xticklabels({'h','J'})

ylim([0 1])


