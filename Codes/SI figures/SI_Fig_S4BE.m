%% ================================ Supplementary Figure 4BE ============================================
close all
clear all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];

% --- Parameters
nbins = 10;

T = [18, 22, 26, 30, 33];
colors = rainbow(numel(T));

% --- Init.
dirisin = dir([isinpath '*.mat']);
is_contra = [];
R_neg = [];

for idf = 1:numel(dirisin)
	
    fnisin = [isinpath dirisin(idf).name];
    D = load(fnisin);
    
    % - DATA
    mLdata = D.Dinference_corr.mleft_data;
    mRdata = D.Dinference_corr.mright_data;
    
    time = linspace(D.Dinference_corr.time(1), D.Dinference_corr.time(end), numel(mLdata));
    loc = strfind(dirisin(idf).name, '.mat');
    idata = str2double(dirisin(idf).name(loc-2:loc-1));
    fn = [dirisin(idf).name(1:3) '_Fish' num2str(idata, '%02i')];
    
    % - ISING
    mLdata = D.Dinference_corr.mleft;
    mRdata = D.Dinference_corr.mright;
    
        
    scat_ising = zeros([size([D.Dinference_corr.leftspikesbin D.Dinference_corr.rightspikesbin]') 3]);
    scat_data = zeros([size([D.Dinference_corr.leftspikesbin_data D.Dinference_corr.rightspikesbin_data]') 3]);
    
    scat_ising(1:length(D.Dinference_corr.L_reg),:,1) = D.Dinference_corr.leftspikesbin';
    scat_ising(length(D.Dinference_corr.L_reg)+1:end,:,3) = D.Dinference_corr.rightspikesbin';
    scat_data(1:length(D.Dinference_corr.L_reg),:,1) = D.Dinference_corr.leftspikesbin_data';
    scat_data(length(D.Dinference_corr.L_reg)+1:end,:,3) = D.Dinference_corr.rightspikesbin_data';
    
    temperatures(idf) = D.Dinference_corr.T;
    T_n = find(temperatures(idf)==T);    
    
    
    R = corrcoef([D.Dinference_corr.leftspikesbin D.Dinference_corr.rightspikesbin]);
    R_data = corrcoef([D.Dinference_corr.leftspikesbin_data D.Dinference_corr.rightspikesbin_data]);
    
    mean_R(T_n,idf,1) = mean(R(D.Dinference_corr.L_reg,D.Dinference_corr.L_reg),[1 2]);
    mean_R(T_n,idf,2) = mean(R(D.Dinference_corr.R_reg,D.Dinference_corr.R_reg),[1 2]);
    mean_R(T_n,idf,3) = mean(R(D.Dinference_corr.R_reg,D.Dinference_corr.L_reg),[1 2]);
    
    mean_R_data(T_n,idf,1) = mean(R_data(D.Dinference_corr.L_reg,D.Dinference_corr.L_reg),[1 2]);
    mean_R_data(T_n,idf,2) = mean(R_data(D.Dinference_corr.R_reg,D.Dinference_corr.R_reg),[1 2]);
    mean_R_data(T_n,idf,3) = mean(R_data(D.Dinference_corr.R_reg,D.Dinference_corr.L_reg),[1 2]);
    
    [ifind jfind] = find(R_data<0);
    Prob_contra_neg(idf) = sum(abs(ismember(ifind,D.Dinference_corr.L_reg)-ismember(jfind,D.Dinference_corr.L_reg))>0)/length(ifind);
    Prob_contra_neg_2(idf) = 0.5*sum(abs(ismember(ifind,D.Dinference_corr.L_reg)-ismember(jfind,D.Dinference_corr.L_reg))>0)/(length(D.Dinference_corr.L_reg)*length(D.Dinference_corr.R_reg));

    thr = 1.2;
    [ifind jfind] = find(R_data<thr);
    is_contra = [is_contra (abs(ismember(ifind,D.Dinference_corr.L_reg)-ismember(jfind,D.Dinference_corr.L_reg))>0)'];
    R_neg = [R_neg R_data(find(R_data<thr))'];
    
    
    
    
    idf
    
end

is_ipsi = abs(1-is_contra);

%%

Y = discretize(R_neg,200);
clear avg_is_contra avg_R_neg
for i = 1:max(Y);
avg_is_contra(i) = mean(is_contra(find(Y==i)));
avg_R_neg(i) = mean(R_neg(find(Y==i)));
end

figure; 
hold on;
plot(avg_R_neg,avg_is_contra,'b','LineWidth',3)


clear avg_is_ipsi avg_R_neg
for i = 1:max(Y);
avg_is_ipsi(i) = mean(is_ipsi(find(Y==i)));
avg_R_neg(i) = mean(R_neg(find(Y==i)));
end

plot(avg_R_neg,avg_is_ipsi,'r','LineWidth',3)
%ylim([0.7 1])
ylabel('Proba to see a contra or an ipsi pair')
xlabel('Correlation')
set(gca,'FontSize',18)

%%


[dist_contra_Y dist_contra_X] = ksdensity(R_neg(is_contra==1));
[dist_ipsi_Y dist_ipsi_X] = ksdensity(R_neg(is_ipsi==1));

figure; hold on
plot(dist_contra_X,dist_contra_Y,'b','LineWidth',2)
plot(dist_ipsi_X,dist_ipsi_Y,'r','LineWidth',2)
xlabel('Correlation')
ylabel('pdf')
set(gca,'FontSize',18)

