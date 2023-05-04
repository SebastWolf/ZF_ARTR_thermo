%% ================================ SI Figure 3G ============================================
%Excess log likelihood ising vs indpdt

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
isinpath_cross_val = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'Cross validation' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

CCC_ALL = [];
CCC_model_ALL = [];
FFF_ALL = [];
FFF_model_ALL = [];
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
        
        id_fish
        t_files
        
        % load data
        

        load([isinpath_cross_val '/T' num2str(T(idT)) '_Fish ' num2str(fish_num(id_fish)) '/ACE_short/_BM_to_hopfield_params_likelihood_MC_save.mat'])
        load([isinpath_cross_val '/T' num2str(T(idT)) '_Fish ' num2str(fish_num(id_fish)) '/ACE_short/_BM_to_hopfield_params_likelihood_on_synth.mat'])
        load([isinpath_cross_val '/T' num2str(T(idT)) '_Fish ' num2str(fish_num(id_fish)) '/ACE_short/_BM_to_hopfield_params_likelihood_on_synth_indpdt.mat'])
        %data_BM = double(data_BM);

        likelihood_on_original_pseudo_ALL(id_fish,t_files) = likelihood_on_original_pseudo;
        likelihood_on_original_save_ALL(id_fish,t_files) = likelihood_on_original_save;
        likelihood_on_test_pseudo_ALL(id_fish,t_files) = likelihood_on_test_pseudo;
        likelihood_on_test_save_ALL(id_fish,t_files) = likelihood_on_test_save;
  
        likelihood_on_test_indpdt_pseudo_ALL(id_fish,t_files) = likelihood_on_synth_pseudo_indpdt;
        likelihood_on_test_indpdt_save_ALL(id_fish,t_files) = likelihood_on_synth_save_indpdt;
        
    end
    
end

%%


D_LL_pseudo = likelihood_on_original_pseudo_ALL-likelihood_on_test_pseudo_ALL;
Eccess_D_LL_pseudo = (likelihood_on_original_pseudo_ALL-likelihood_on_test_pseudo_ALL)./likelihood_on_test_pseudo_ALL;
D_LL_me = likelihood_on_original_save_ALL-likelihood_on_test_save_ALL;

D_LL_me(D_LL_me==0) = NaN;
D_LL_pseudo(D_LL_pseudo==0) = NaN;
Eccess_D_LL_pseudo

%%

D_LL_pseudo = likelihood_on_original_pseudo_ALL-likelihood_on_test_pseudo_ALL;
D_LL_me = likelihood_on_original_save_ALL-likelihood_on_test_save_ALL;

Eccess_D_LL_pseudo = (likelihood_on_original_pseudo_ALL-likelihood_on_test_indpdt_pseudo_ALL)./likelihood_on_test_indpdt_pseudo_ALL;
Eccess_D_LL_test_pseudo = (likelihood_on_test_pseudo_ALL-likelihood_on_test_indpdt_pseudo_ALL)./likelihood_on_test_indpdt_pseudo_ALL;

D_LL_me(D_LL_me==0) = NaN;
D_LL_pseudo(D_LL_pseudo==0) = NaN;
Eccess_D_LL_pseudo


%%

figure;
boxplot(-[1*(Eccess_D_LL_pseudo(:)) 1*(Eccess_D_LL_test_pseudo(:))])
set(gca,'FontSize',18)
box off
xlabel('')

ylim([0 1])
