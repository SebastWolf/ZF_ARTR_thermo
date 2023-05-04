%% ================================ SI Figure 6 inset ============================================

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

for id_fish = 4
    clear J_ h_
    dpath = hpath_star(fish_num(id_fish));
    clear files_T
    files = dir(dpath);
    count_1 = 0;
    
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
        T_(t_files) = T(idT);
        
    end
    
    
    for t_files1 = 1:length(files)
        count_1 = count_1+1;
        for t_files2 = t_files1+1:length(files)
            
            X_J = squeeze(J_(:,:,t_files1));
            Y_J = squeeze(J_(:,:,t_files2));
            
            X_h = squeeze(h_(:,t_files1));
            Y_h = squeeze(h_(:,t_files2));
            
            figure
            plot(X_J,Y_J,'k.','MarkerSize',10)
            title([files(t_files).name ' - TX1 = ' num2str(T_(t_files1)) ' - TX2 = ' num2str(T_(t_files2))])
            xlabel('J TX1')
            ylabel('J TX2')
            set(gca,'FontSize',18)
            %set(gca, 'YScale', 'log')
            %set(gca, 'XScale', 'log')
            xlim([-1, 2]); ylim([-1, 2]);
            box off

            
            figure
            plot(X_h,Y_h,'k.','MarkerSize',10)
            title([files(t_files).name ' - TX1 = ' num2str(T_(t_files1)) ' - TX2 = ' num2str(T_(t_files2))])
            xlabel('h TX1')
            ylabel('h TX2')
            set(gca,'FontSize',18)
            xlim([-8, 0]); ylim([-8, 0]);
            box off

            
            
            R = corrcoef([X_J(:),Y_J(:)]);
            R_J_keep = R(1,2);
            
            R = corrcoef([X_h(:),Y_h(:)]);
            R_h_keep = R(1,2);
            
            
            figure;
            h = histogram2(X_J,Y_J,'DisplayStyle','tile','ShowEmptyBins','off');
            colorbar
            title([files(t_files).name ' - TX1 = ' num2str(T_(t_files1)) ' - TX2 = ' num2str(T_(t_files2)) ' R2 = ' num2str(R_J_keep^2)])
            xlim([-1, 2]); ylim([-1, 2]);
            xlabel('J TX1')
            ylabel('J TX2')
            set(gca,'FontSize',18)
            box off
            
            figure;
            h = histogram2(X_h,Y_h,'DisplayStyle','tile','ShowEmptyBins','off');
            colorbar
            title([files(t_files).name ' - TX1 = ' num2str(T_(t_files1)) ' - TX2 = ' num2str(T_(t_files2)) ' R2 = ' num2str(R_h_keep^2)])
            xlabel('h TX1')
            ylabel('h TX2')
            set(gca,'FontSize',18)
            xlim([-10, 2]); ylim([-8, 0]);
            box off
 
            
        end
    end
    
end


%%


figure;
boxplot([R_h_keep.^2' R_J_keep.^2'])
set(gca,'FontSize',18)
box off
xticks([1 2])
xticklabels({'h','J'})

ylim([0 1])

