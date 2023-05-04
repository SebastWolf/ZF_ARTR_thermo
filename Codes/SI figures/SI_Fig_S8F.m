%% ================================ Supplementary Figure 8F ============================================

clear all
close all
clc

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
hpath_light = [pwd filesep 'Data availability' filesep 'Light stim' filesep ];

dpath = hpath_light;
files = dir([dpath '/*.mat']);
nfiles = numel(files);

list_idf = [1:numel(files)]
pcount = 0;

int_spont_LL_ = [];
int_spont_RR_ = [];
int_spont_LR_ = [];
int_R_LR_ = [];
int_L_LR_ = [];
int_R_RR_ = [];
int_R_LL_ = [];
int_L_LL_ = [];
int_L_RR_ = [];


for idf = list_idf
    idf
    load([files(idf).folder filesep files(idf).name])   

    pcount = pcount+1;
    A_cor_rh23_ALL = Dinference.A_cor_rh23_ALL;
    L_vec = Dinference.L_reg;
    R_vec = Dinference.R_reg;
    CONFIG_h_sup = Dinference.CONFIG_under_stim_LONG;
    vect_stim = Dinference.vect_stim_LONG;
    
    
    time_stim_L = find(vect_stim==1);
    time_stim_R = find(vect_stim==-1);
    
    time_of_interest = time_stim_R;
    AA = CONFIG_h_sup(time_of_interest,:);
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_R_LL_all(idf) = nanmean(int(:)); int_R_LL_ = [int_R_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_R_RR_all(idf) = nanmean(int(:)); int_R_RR_ = [int_R_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_R_LR_all(idf) = nanmean(int(:)); int_R_LR_ = [int_R_LR_  int(:)'];
    
    time_of_interest = time_stim_L;
    AA = CONFIG_h_sup(time_of_interest,:);
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_L_LL_all(idf) = nanmean(int(:)); int_L_LL_ = [int_L_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_L_RR_all(idf) = nanmean(int(:)); int_L_RR_ = [int_L_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_L_LR_all(idf) = nanmean(int(:)); int_L_LR_ = [int_L_LR_  int(:)'];
    
    AA = Dinference.CONFIG;
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_spont_LL_all(idf) = nanmean(int(:)); int_spont_LL_ = [int_spont_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_spont_RR_all(idf) = nanmean(int(:)); int_spont_RR_ = [int_spont_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_spont_LR_all(idf) = nanmean(int(:)); int_spont_LR_ = [int_spont_LR_  int(:)'];
   
end

int_spont_LL_(int_spont_LL_==1) = NaN;
int_spont_LR_(int_spont_LR_==1) = NaN;
int_R_LL_(int_R_LL_==1) = NaN;
int_R_RR_(int_R_RR_==1) = NaN;
int_L_LL_(int_L_LL_==1) = NaN;
int_L_RR_(int_L_RR_==1) = NaN;
int_R_LR_(int_R_LR_==1) = NaN;
int_L_LR_(int_L_LR_==1) = NaN;

C_spont_LL_all(C_spont_LL_all==0) = NaN;
C_spont_RR_all(C_spont_RR_all==0) = NaN;
C_spont_LR_all(C_spont_LR_all==0) = NaN;

C_L_RR_all(C_L_RR_all==0) = NaN;
C_L_LL_all(C_L_LL_all==0) = NaN;

C_R_LL_all(C_R_LL_all==0) = NaN;
C_R_RR_all(C_R_RR_all==0) = NaN;
C_L_LL_all(C_L_LL_all==0) = NaN;

C_L_LL_all(C_L_LL_all==0) = NaN;
C_L_RR_all(C_L_RR_all==0) = NaN;
C_L_LR_all(C_L_LR_all==0) = NaN;
C_R_LR_all(C_R_LR_all==0) = NaN;


%% ================================ Supplementary Figure 8F ============================================

figure;
boxplot([(C_R_LR_all+C_L_LR_all)'/2 C_spont_LR_all'])
set(gca,'FontSize',18)
box off
xticks([1:3])
xticklabels({'Under stimulation','Spont'})
ylabel('<R_{i,j}>_{contra}')

figure;
boxplot([(C_L_LL_all+C_R_RR_all+C_R_LL_all+C_L_RR_all)'/4  (C_spont_RR_all+C_spont_LL_all)'/2])
ylim([-0 0.6 ])
set(gca,'FontSize',18)
box off
xticks([1:2])
xticklabels({'Under stimulation','Under spont'})
ylabel('<R_{i,j}>_{ipsi}')


