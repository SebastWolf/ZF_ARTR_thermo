%% ================================ Figure 5D ============================================
% Values of the additional biases averaged over the ipsilateral and contralateral (with respect to the stimulated eye) neural populations.

clear all
close all
clc

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
hpath_light = [pwd filesep 'Data availability' filesep 'Light stim' filesep ];

files = dir([hpath_light '/*.mat']);
pcount = 0;

list_idf = [1:numel(files)]

for idf = list_idf
    idf
    load([files(idf).folder filesep files(idf).name])   

    pcount = pcount+1;
    
    L_vec = Dinference.L_reg;
    R_vec = Dinference.R_reg;    

    h_sup_R = Dinference.h_sup_R;
    h_sup_L = Dinference.h_sup_L;
    h_sup_spont = Dinference.h_sup_spont;

    
    H_RL_all(idf) = nanmean(h_sup_R(L_vec))
    H_RR_all(idf) =  nanmean(h_sup_R(R_vec))
    
    H_LL_all(idf) = nanmean(h_sup_L(L_vec))
    H_LR_all(idf) = nanmean(h_sup_L(R_vec))
    
    H_Lspont_all(idf) = nanmean(h_sup_spont(L_vec))
    H_Rspont_all(idf) = nanmean(h_sup_spont(R_vec))
    
end

%%

H_RR_all(H_RR_all==0) = NaN;
H_LR_all(H_LR_all==0) = NaN;
H_LL_all(H_LL_all==0) = NaN;
H_RL_all(H_RL_all==0) = NaN;

figure;
boxplot([(H_RR_all+H_LL_all)'/2 (H_RL_all+H_LR_all)'/2])
ylim([-0.1 0.5 ])
set(gca,'FontSize',16)
box off
xticks([1 2])

xticklabels({'Ipsi','Contra'})
ylabel('<H>')
