 %% =========================== Supplementary Figure 7E ============================================

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

pcount = 0;

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
        
        J_J_L = table_mean_field(line_selected,3);
        J_J_R = table_mean_field(line_selected,4);
        I_I = table_mean_field(line_selected,5);
        M_h_L = table_mean_field(line_selected,6);
        M_h_R = table_mean_field(line_selected,7);
        K_L_all_fit(idT,id_fish) =  table_mean_field(line_selected,8);
        K_R_all_fit(idT,id_fish) =  table_mean_field(line_selected,9);
        
        pcount = pcount + 1;
        J_L_all_data(idT,id_fish) = J_J_L;
        J_R_all_data(idT,id_fish) = J_J_R;
        I_all_data(idT,id_fish) = I_I;
        M_h_R_all_data(idT,id_fish) = M_h_R;
        M_h_L_all_data(idT,id_fish) = M_h_L;
        
        fish_num_all_data(idT,id_fish) = fish_num(id_fish);
        fish_temp_all_data(idT,id_fish) = T(idT);
       
    end 
end

%%

J_R_all_data(J_R_all_data==0) = NaN;
J_L_all_data(J_L_all_data==0) = NaN;
M_h_R_all_data(M_h_R_all_data==0) = NaN;
M_h_L_all_data(M_h_L_all_data==0) = NaN;
I_all_data(I_all_data==0) = NaN;

cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

for idT = 1:numel(T)
figure(1000); hold on
plot(M_h_R_all_data(idT,:),J_R_all_data(idT,:),'.','MarkerSize',20, 'MarkerEdgeColor',cmp(idT,:))
plot(M_h_L_all_data(idT,:),J_L_all_data(idT,:),'.','MarkerSize',20, 'MarkerEdgeColor',cmp(idT,:))
end
figure(1000);
set(gca,'FontSize',18)
xlabel('H')
ylabel('J')

%%

M_H_LIST = [M_h_R_all_data(:)' M_h_L_all_data(:)']
J_LIST = [J_R_all_data(:)' J_L_all_data(:)']
M_H_LIST(isnan(M_H_LIST)==1)=[];
J_LIST(isnan(J_LIST)==1)=[];

fir = fit(M_H_LIST',J_LIST','poly1')
p1 = fir.p1;
p2 = fir.p2;
plot(M_H_LIST,fir(M_H_LIST))

%%

J_R_all_data(J_R_all_data==0) = NaN;
J_L_all_data(J_L_all_data==0) = NaN;
M_h_R_all_data(M_h_R_all_data==0) = NaN;
M_h_L_all_data(M_h_L_all_data==0) = NaN;
I_all_data(I_all_data==0) = NaN;

cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

for idT = 1:numel(T)
figure(1000); hold on
plot(M_h_R_all_data(idT,:),J_R_all_data(idT,:),'.','MarkerSize',20, 'MarkerEdgeColor',cmp(idT,:))
plot(M_h_L_all_data(idT,:),J_L_all_data(idT,:),'.','MarkerSize',20, 'MarkerEdgeColor',cmp(idT,:))
end
figure(1000);
set(gca,'FontSize',18)
xlabel('H')
ylabel('J')

plot([-5:0.1:-3.5],fir([-5:0.1:-3.5]),'r','LineWidth',3)

%%

Ka = -1/fir.p1;

clear var_L var_R
for id_fish = 1:numel(fish_num)
    var_L(:,id_fish) = K_L_all_fit(:,id_fish).*(M_h_L_all_data(:,id_fish)+Ka*J_L_all_data(:,id_fish));
    var_R(:,id_fish) = K_R_all_fit(:,id_fish).*(M_h_R_all_data(:,id_fish)+Ka*J_R_all_data(:,id_fish));
end


%%
var_all = cat(2,var_L,var_R);
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
