%% ================================ Figure 5FG ============================================
% Mean-field energy landscape with and without additional inputs

clear all
close all
clc

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
hpath_light = [pwd filesep 'Data availability' filesep 'Light stim' filesep ];
load([pwd '/Data availability/Mean_field_light_Parameters.mat'])

cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

dpath = hpath_light;
files = dir([dpath '/*.mat']);
nfiles = numel(files);
pcount = 0;
list_idf = [1:numel(files)]

for i=1:6
idf = list_idf(i)
line_selected = idf
idf
load([files(idf).folder filesep files(idf).name])

pcount = pcount+1;


J_MC = Dinference.J_MC;
h_MC = Dinference.h_MC;
J = J_MC;
h = h_MC';

L_reg = Dinference.L_reg;
R_reg = Dinference.R_reg;
N_L = length(L_reg);
N_R = length(R_reg);
J_inter = Dinference.J_MC(1:N_L,N_L+1:N_L+N_R);
J_left = Dinference.J_MC(1:N_L,1:N_L);
J_right = Dinference.J_MC(N_L+1:N_L+N_R,N_L+1:N_L+N_R) ;


%% ================== . Free energy without input


M_h_R = table_mean_field_light(line_selected,7);
M_h_L = table_mean_field_light(line_selected,6);;

M_J_R = nanmean(J_right(:))
M_J_L = nanmean(J_left(:))
M_I_I = nanmean(J_inter(:));

J_J_L = table_mean_field_light(line_selected,3);
J_J_R = table_mean_field_light(line_selected,4);
I_I = table_mean_field_light(line_selected,5);

% activity

N_L_noise = table_mean_field_light(line_selected,8);
N_R_noise = table_mean_field_light(line_selected,9);

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

figure
sss = surf(free_energy_all)
sss.EdgeColor = 'none';
title(['SPONT'])
set(gca,'FontSize',18)


% ENERGY PATH
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
pcount = pcount+1;
N_points_0 = find(diff(stat_points')<=length(u_list)/2 & sum(stat_points')<=length(u_list)/2);
if isempty(find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2))==0;
    N_points_1_2 = find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2);
    DE_0_12(idf) = energy_stat_points(N_points_1_2)-energy_stat_points(N_points_0);
    
    [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_1_2,:), stat_points(N_points_0,:), free_energy_all);
    [path_X path_Y] = ind2sub(size(free_energy_all),path);
    
    clear energy_along
    for ttt = 1:length(path_X)
        energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
    end
    figure(idf+30); hold on;
    plot((path_Y-path_X)/length(u_list),energy_along,'Color',cmp(1,:))
    
    Delta_0_12(idf,1) = max(energy_along)-energy_along(end); %0 vers 1
    Delta_0_12(idf,2) = max(energy_along)-energy_along(1); %1 vers 0
    
end



if isempty(find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2))==0;
    N_points_2_1 = find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2);
    DE_0_21(idf) = energy_stat_points(N_points_2_1)-energy_stat_points(N_points_0);
    
    [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_2_1,:), stat_points(N_points_0,:), free_energy_all);
    [path_X path_Y] = ind2sub(size(free_energy_all),path);
    
    clear energy_along
    for ttt = 1:length(path_X)
        energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
    end
    figure(idf+30); hold on;
    plot((path_Y-path_X)/length(u_list),energy_along,'Color',cmp(1,:))
    set(gca,'FontSize',16)
    
    Delta_0_21(idf,1) = max(energy_along)-energy_along(end); %0 vers 1
    Delta_0_21(idf,2) = max(energy_along)-energy_along(1); %1 vers 0
    
end
title(['SPONT'])
set(gca,'FontSize',18)



%% ================ Free energy with input


H_LL_all(idf) = nanmean(Dinference.h_sup_L(L_reg))
H_LR_all(idf) = nanmean(Dinference.h_sup_L(R_reg))


M_h_R = table_mean_field_light(line_selected,7)+H_LR_all(idf);
M_h_L = table_mean_field_light(line_selected,6)+H_LL_all(idf);

M_J_R = nanmean(J_right(:))
M_J_L = nanmean(J_left(:))
M_I_I = nanmean(J_inter(:));

J_J_L = table_mean_field_light(line_selected,3);
J_J_R = table_mean_field_light(line_selected,4);
I_I = table_mean_field_light(line_selected,5);

% activity

N_L_noise = table_mean_field_light(line_selected,8);
N_R_noise = table_mean_field_light(line_selected,9);

N_L = N_L_noise;
N_R = N_R_noise;
u_list = [1e-4:0.01:1-1e-4];
v_list = [1e-4:0.01:1-1e-4];


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
title(['WITH INPUT'])
set(gca,'FontSize',18)


% ENERGY PATH
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
pcount = pcount+1;
N_points_0 = find(diff(stat_points')<=length(u_list)/2 & sum(stat_points')<=length(u_list)/2);
if isempty(find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2))==0;
    N_points_1_2 = find(stat_points(:,1)>=1*length(u_list)/2 & stat_points(:,2)<=1*length(u_list)/2);
    DE_0_12(idf) = energy_stat_points(N_points_1_2)-energy_stat_points(N_points_0);
    
    [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_1_2,:), stat_points(N_points_0,:), free_energy_all);
    [path_X path_Y] = ind2sub(size(free_energy_all),path);
    
    clear energy_along
    for ttt = 1:length(path_X)
        energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
    end
    figure(idf+30); hold on;
    plot((path_Y-path_X)/length(u_list),energy_along,'Color',cmp(3,:))
    
    Delta_0_12(idf,1) = max(energy_along)-energy_along(end); %0 vers 1
    Delta_0_12(idf,2) = max(energy_along)-energy_along(1); %1 vers 0
    
end

if isempty(find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2))==0;
    N_points_2_1 = find(stat_points(:,1)<=1*length(u_list)/2 & stat_points(:,2)>=1*length(u_list)/2);
    DE_0_21(idf) = energy_stat_points(N_points_2_1)-energy_stat_points(N_points_0);
    
    [dist,path,pred] = shortest_path_free_energy(stat_points(N_points_2_1,:), stat_points(N_points_0,:), free_energy_all);
    [path_X path_Y] = ind2sub(size(free_energy_all),path);
    
    clear energy_along
    for ttt = 1:length(path_X)
        energy_along(ttt) = free_energy_all(path_X(ttt),path_Y(ttt));
    end
    figure(idf+30); hold on;
    plot((path_Y-path_X)/length(u_list),energy_along,'Color',cmp(3,:))
    set(gca,'FontSize',16)
    
    Delta_0_21(idf,1) = max(energy_along)-energy_along(end); %0 vers 1
    Delta_0_21(idf,2) = max(energy_along)-energy_along(1); %1 vers 0
    
end

title(['WITH INPUT'])
set(gca,'FontSize',18)

end