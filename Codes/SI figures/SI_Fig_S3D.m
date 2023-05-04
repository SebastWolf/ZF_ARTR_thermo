%% =========================== Supplementary Figure 3D ============================================


% Probability that K of the N neurons in the ARTR are active simultaneously in the data (black dots) and in the model (yellow line).
% To obtain these probabilities, we pool together all states where any K neurons are active together while the rest of the
% neurons are silent.

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

id_fish = 6

dpath = hpath_star(fish_num(id_fish));
files = dir(dpath);
t_files = 1;
clear files_T

for idT = 1:numel(T)
    if isempty(strfind(files(t_files).name,num2str(T(idT))))==0;
        files_T(t_files) = T(idT);
    end
end
idT = find(files_T(t_files)==T);

id_fish
t_files

load([files(t_files).folder filesep files(t_files).name]);

Nneurons = length(Dinference_corr.h);

data_expe = [Dinference_corr.leftspikesbin_data Dinference_corr.rightspikesbin_data];
data_BM = [Dinference_corr.leftspikesbin Dinference_corr.rightspikesbin];
clear pk
config = data_expe';
number_of_active = sum(config,1);
for k = 1:Nneurons;
    pk(k) = length(find(number_of_active==k))/length(number_of_active);
end

clear pk_sim
config = data_BM;
number_of_active = sum(config,2);
for k = 1:Nneurons;
    pk_sim(k) = length(find(number_of_active==k))/length(number_of_active);
end

figure('Name', files(t_files).name); hold on
hold on;
plot(pk_sim,'r','LineWidth',3)
plot(pk,'ko','LineWidth',3)
set(gca,'YScale','log')
set(gca,'FontSize',18)
xlim([0 225])
