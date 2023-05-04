%% ================================ Supplementary Figure 8CDE ============================================

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

fHz = 1;
Tfinal = 3990;
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
    config = A_cor_rh23_ALL;

    %stimulation times
    
    DeltaT=[30 10 25 5 20 15];
    Tinit=[10 1100 1530 2455 2720 3480];
    Tfinal=[1000 1440 2380 2625 3400 3990];
    
    time_spont = [];
    for p = 1:length(DeltaT)-1;
        time_spont = [time_spont Tfinal(p)+DeltaT(p)+1:Tinit(p+1)-1];
    end
    time_spont = [time_spont Tfinal(length(Tfinal)):size(A_cor_rh23_ALL,1)];
        
    out_delt = find(DeltaT==5);
    DeltaT(out_delt) = []; Tinit(out_delt) = []; Tfinal(out_delt) = [];
    
    dec = 5;    
    ind_flash_L = [];
    ind_flash_R = [];    
    
    for p = 1:length(DeltaT)
        
        ind_flash_R_int = [Tinit(p)*fHz+DeltaT(p)*fHz+2:2*DeltaT(p)*fHz:Tfinal(p)*fHz];
        ind_flash_L_int = [Tinit(p)*fHz+2:2*DeltaT(p)*fHz:Tfinal(p)*fHz];
        %
       
        for lll = 1:length(ind_flash_R_int);
            ind_flash_R = [ind_flash_R ind_flash_R_int(lll):ind_flash_R_int(lll)+DeltaT(p)*fHz];
        end
        
        for lll = 1:length(ind_flash_L_int);
            ind_flash_L = [ind_flash_L ind_flash_L_int(lll):ind_flash_L_int(lll)+DeltaT(p)*fHz];
        end           
    end   
    
    time_stim_L = ind_flash_L;
    time_stim_R = ind_flash_R;
    
    time_stim_L = unique(time_stim_L);
    time_stim_R = unique(time_stim_R);
    
    vect_stim_data_L = 0*config(:,1);
    vect_stim_data_L(time_stim_L) = 1;
    
    vect_stim_data_R = 0*config(:,1);
    vect_stim_data_R(time_stim_R) = 1;
   
    time_of_interest = time_stim_R;
    AA = config(time_of_interest,:);
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_R_LL_all(idf) = nanmean(int(:)); int_R_LL_ = [int_R_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_R_RR_all(idf) = nanmean(int(:)); int_R_RR_ = [int_R_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_R_LR_all(idf) = nanmean(int(:)); int_R_LR_ = [int_R_LR_  int(:)'];
    
    time_of_interest = time_stim_L;
    AA = config(time_of_interest,:);
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_L_LL_all(idf) = nanmean(int(:)); int_L_LL_ = [int_L_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_L_RR_all(idf) = nanmean(int(:)); int_L_RR_ = [int_L_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_L_LR_all(idf) = nanmean(int(:)); int_L_LR_ = [int_L_LR_  int(:)'];
    
    time_of_interest = time_spont;
    AA = config(time_of_interest,:);
    CC_AA = corrcoef(AA);
    int = CC_AA(L_vec,L_vec); C_spont_LL_all(idf) = nanmean(int(:)); int_spont_LL_ = [int_spont_LL_  int(:)'];
    int = CC_AA(R_vec,R_vec); C_spont_RR_all(idf) = nanmean(int(:)); int_spont_RR_ = [int_spont_RR_  int(:)'];
    int = CC_AA(L_vec,R_vec); C_spont_LR_all(idf) = nanmean(int(:)); int_spont_LR_ = [int_spont_LR_  int(:)'];
    
    close all
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
C_L_LR_all(C_L_LR_all==0) = NaN;


%% ================================ Supplementary Figure 8C ============================================
xmin = -1;
xmax = 1;

figure;
plot(int_spont_LR_(:),int_R_LR_(:),'k.')
hold on
plot(int_spont_LR_(:),int_L_LR_(:),'k.')

plot([xmin xmax],[xmin xmax],'r')
set(gca,'FontSize',18)
xlabel('spont')
ylabel('stim')
box off
xlim([xmin xmax])
ylim([xmin xmax])

%% ================================ Supplementary Figure 8D ============================================

figure;
boxplot([(C_R_LR_all+C_L_LR_all)'/2 C_spont_LR_all'])
set(gca,'FontSize',18)
box off
xticks([1:2])
xticklabels({'Under stimulation','Spont'})
ylabel('<R_{i,j}>_{contra}')
ylim([-0.2 0.6 ])

%%
ct1 = (C_R_LR_all+C_L_LR_all)'/2;
ct2 = C_spont_LR_all';

[h,p] = ttest(ct1,ct2)
%%
figure;
boxplot([(C_L_LL_all+C_R_RR_all+C_R_LL_all+C_L_RR_all)'/4  (C_spont_RR_all+C_spont_LL_all)'/2])
set(gca,'FontSize',18)
box off
xticks([1:2])
xticklabels({'Under stimulation','Under spont'})
ylabel('<R_{i,j}>_{ipsi}')
ylim([-0.2 0.6 ])


ct1 = (C_L_LL_all+C_R_RR_all+C_R_LL_all+C_L_RR_all)'/4;
ct2 = (C_spont_RR_all+C_spont_LL_all)'/2;

[h,p] = ttest(ct1,ct2)

%% ================================ Supplementary Figure 8E ============================================

xmin = -0.5;
xmax = 1;

figure; hold on
plot(int_spont_LL_(:),int_L_LL_(:),'k.')
plot(int_spont_LL_(:),int_R_LL_(:),'k.')
plot(int_spont_RR_(:),int_R_RR_(:),'k.')
plot(int_spont_RR_(:),int_L_RR_(:),'k.')
hold on
plot([xmin xmax],[xmin xmax],'r')
set(gca,'FontSize',18)
xlabel('spont')
ylabel('stim')
box off
xlim([xmin xmax])
ylim([xmin xmax])
