%% ================================ Supplementary Figure 8AB ============================================
% To assess the quality of the model of the ARTR of these visually-driven experiments, we compare the mean 
% activity (C) and the pairwise covariance (D) computed on real data (spontaneous part of the recordings) to
% those computed on the synthetic data.

clear all
close all
clc

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
hpath_light = [pwd filesep 'Data availability' filesep 'Light stim' filesep ];

dpath = hpath_light;
files = dir([dpath '/*.mat']);
nfiles = numel(files);
pcount = 0;
CCC_ALL = [];
CCC_model_MH_ALL = [];
FFF_ALL = [];
FFF_model_MH_ALL = [];
list_idf = [1:numel(files)]

for idf = list_idf
    idf
    load([files(idf).folder filesep files(idf).name])   

    pcount = pcount+1;
    A_cor_rh23 = Dinference.A_cor_rh23;
    L_vec = Dinference.L_reg;
    R_vec = Dinference.R_reg;
    vect_stim = Dinference.vect_stim_LONG;
    CONFIG = Dinference.CONFIG;
   
    CCC = cov(A_cor_rh23);
    CCC_ALL = [CCC_ALL CCC(:)'];
    
    FFF = nanmean(A_cor_rh23,1);
    FFF_ALL = [FFF_ALL FFF];
    
    CCC_model_MH = cov(CONFIG);
    CCC_model_MH_ALL = [CCC_model_MH_ALL CCC_model_MH(:)'];

    FFF_model_MH = nanmean(CONFIG,1);
    FFF_model_MH_ALL = [FFF_model_MH_ALL FFF_model_MH];
    
end

fit_FFF_MH = fit(FFF_ALL',FFF_model_MH_ALL','poly1');
fit_CCC_MH = fit(CCC_ALL',CCC_model_MH_ALL','poly1');

%% ================================ Supplementary Figure 8A ============================================

xmin = 0;
xmax = 1;
figure; clf
hold on;
plot(FFF_ALL,FFF_model_MH_ALL,'k.')
plot([0 1],fit_FFF_MH([0 1]),'r','LineWidth',1)
xlabel('data')
ylabel('model')
xlim([xmin xmax])
ylim([xmin xmax])
set(gca,'FontSize',18)

%% ================================ Supplementary Figure 8B ============================================

xmin = -0.3;
xmax = 0.4;
figure; clf
hold on;
plot(CCC_ALL,CCC_model_MH_ALL,'k.')
plot([xmin xmax]',fit_CCC_MH([xmin xmax]'),'r','LineWidth',1)
xlabel('data')
ylabel('model')
set(gca,'FontSize',18)
xlim([xmin xmax])
ylim([xmin xmax])


