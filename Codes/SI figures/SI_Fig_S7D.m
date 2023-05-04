%% =========================== Supplementary Figure 7D ============================================
% Free-energy difference between stationary sates of the landscape as a function of the temperature.

clc
clear all
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];
load([pwd filesep 'Data availability' filesep 'Mean_field_Parameters.mat'])

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

for id_fish = 1:numel(fish_num)
    
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
        
        %% compute landscape
       
        M_h_R = table_mean_field(line_selected,7);
        M_h_L = table_mean_field(line_selected,6);;

        M_J_R = nanmean(J_right(:))       
        M_J_L = nanmean(J_left(:))        
        M_I_I = nanmean(J_inter(:));
        
        J_J_L = table_mean_field(line_selected,3);
        J_J_R = table_mean_field(line_selected,4);
        I_I = table_mean_field(line_selected,5);
        
        % activity
                        
        N_L_noise = table_mean_field(line_selected,8);
        N_R_noise = table_mean_field(line_selected,9);
                
        N_L = N_L_noise;
        N_R = N_R_noise;
        
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
          
         BW = imregionalmin(free_energy_all);
        [xm ym] = find(imregionalmin(free_energy_all)==1);
        
        clear var
        for i = 1:length(free_energy_all(:,1)); var(i) = free_energy_all(i,i); end
        nem = diff(var); tes = find(nem(2:end).*nem(1:end-1)<0 & diff(diff(var))>0);
        av = (xm+ym)/2;
        for l_tes = 1:length(tes);
            if ismember(tes(l_tes),xm)==0 & max(abs(tes(l_tes)-av)>length(u_list)*0.6)==1 & min(abs(tes(l_tes)-av))>length(u_list)*0.5; xm = cat(1,xm,tes(l_tes)); ym = cat(1,ym,tes(l_tes)); end
        end
        stat_points = [xm ym];
        N_stat_points = length(xm);
        
        
        %extract deltaE
        clear frrrr;
        for N_points = 1:length(xm); frrrr(N_points) = free_energy_all(xm(N_points),ym(N_points)); end
        energy_stat_points = frrrr ;

        N_points_0 = find(diff(stat_points')<=length(u_list)/2 & sum(stat_points')<=length(u_list)/2);
        if isempty(find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2))==0;
            N_points_1_2 = find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2);
            DE_0_12(idT,id_fish) = energy_stat_points(N_points_1_2)-energy_stat_points(N_points_0);
            
            [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_1_2,:), stat_points(N_points_0,:), free_energy_all);
            [path_X path_Y] = ind2sub(size(free_energy_all),path);
            
            clear energy_along
            for ttt = 1:length(path_X)
                energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
            end
             
            Delta_0_12(idT,id_fish,1) = max(energy_along)-energy_along(end); %0 vers 1
            Delta_0_12(idT,id_fish,2) = max(energy_along)-energy_along(1); %1 vers 0
            
        end
        
        
               
        if isempty(find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2))==0;
            N_points_2_1 = find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2);
            DE_0_21(idT,id_fish) = energy_stat_points(N_points_2_1)-energy_stat_points(N_points_0);
            
            [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_2_1,:), stat_points(N_points_0,:), free_energy_all);
            [path_X path_Y] = ind2sub(size(free_energy_all),path);
            
            clear energy_along
            for ttt = 1:length(path_X)
                energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
            end
           
            Delta_0_21(idT,id_fish,1) = max(energy_along)-energy_along(end); %0 vers 1
            Delta_0_21(idT,id_fish,2) = max(energy_along)-energy_along(1); %1 vers 0
            
        end
        
  close all
  
    end
    
end

%%
DE_0_12(DE_0_12==0) = NaN;
DE_0_21(DE_0_21==0) = NaN;
Delta_0_12(Delta_0_12==0) = NaN
Delta_0_12(Delta_0_21==0) = NaN



var_all = cat(2,DE_0_21,DE_0_12);
Y_var = nanmean(var_all,2);
Y_std = nanstd(var_all');

figure; hold on
errorbar(Y_var,nanstd(var_all')./sqrt(sum(isnan(var_all)==0,2))','k','LineWidth',2)

for idT = 1:numel(T)
plot(idT,Y_var(idT),'.','MarkerSize',30, 'MarkerEdgeColor',cmp(idT,:))
end

xticks([1:5])
xticklabels({'18','22','26','30','33'})
set(gca,'FontSize',16)
box off
xlabel('T (°C)')
xlim([1,5])



