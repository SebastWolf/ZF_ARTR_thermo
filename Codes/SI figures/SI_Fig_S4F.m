%% ================================ Supplementary Figure 4F ============================================

clear all
clc
close all


pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];

T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

dist_HBO_ALL = [];
J_ALL = [];
J_L_ALL = [];
J_R_ALL = [];
J_LR_ALL = [];

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
        
        % load Dinference_corr
        load([files(t_files).folder filesep files(t_files).name]);
        
        h_ = Dinference_corr.h;
        J_ = Dinference_corr.J;
        J_ALL = [J_ALL J_(:)'];
        dist_HBO = Dinference_corr.dist_HBO;
        dist_HBO_ALL = [dist_HBO_ALL dist_HBO(:)'];

        
    end
end

%%
dis = dist_HBO_ALL;
conn = J_ALL;
dis(conn==0) = [];
conn(conn==0) = [];

pas = 0.01;
vet = [min(dis(:)):pas:max(dis(:))];
clear conn_moyenne conn_std
for i = 1:length(vet)-1;
    conn_moyenne(i) = mean(conn(find(dis(:)>vet(i) & dis(:)<=vet(i+1))));
    conn_std(i) = std(conn(find(dis(:)>vet(i) & dis(:)<=vet(i+1))));%/sqrt(length((find(dis(:)>vet(i) & dis(:)<=vet(i+1)))));
end


figure(100); clf
options.handle     = figure(100);
options.color_area = [.8 .8 .5]/255;
options.color_line = [1 0 0]/255;
options.alpha = 0.2;
options.line_width = 2; 
options.x_axis = 1e3*(vet(1:end-1)+0.5*pas)  ;    
options.error = 'std' 

data_mean = conn_moyenne;
data_std = conn_std;
plot_areaerrorbar_mod(data_mean, data_std, options)
box off
xlabel('Dij (in um)')
ylabel('Connectivity')
set(gca,'FontSize',18)

