%% ================================ SI Figure 2C ============================================
% ------- average correlation vs temperature

clear all
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
fnout = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep '/Behavior data/meanmLmR.mat'];

% --- Parameters
nbins = 10;

T = [18, 22, 26, 30, 33];
colors = rainbow(numel(T));

% --- Init.
dirisin = dir([isinpath '*.mat']);

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
    mean_R_data(T_n,idf,4) = 0.5*(mean(R_data(D.Dinference_corr.R_reg,D.Dinference_corr.R_reg),[1 2])+mean(R_data(D.Dinference_corr.L_reg,D.Dinference_corr.L_reg),[1 2]));
    
    
    idf
    
end

mean_R(mean_R==0) = NaN;
mean_R_data(mean_R_data==0) = NaN;

%%

figure; hold on;
errorbar(T,nanmean(squeeze(mean_R_data(:,:,1)),2),nanstd(squeeze(mean_R_data(:,:,1))'),'r*')
errorbar(T,nanmean(squeeze(mean_R_data(:,:,2)),2),nanstd(squeeze(mean_R_data(:,:,2))'),'b*')
errorbar(T,nanmean(squeeze(mean_R_data(:,:,3)),2),nanstd(squeeze(mean_R_data(:,:,3))'),'k*')
plot(T,nanmean(squeeze(mean_R_data(:,:,1)),2),'r')
plot(T,nanmean(squeeze(mean_R_data(:,:,2)),2),'b')
plot(T,nanmean(squeeze(mean_R_data(:,:,3)),2),'k')
title('Correlation on data')
ylim([-0.2 0.5])
set(gca,'FontSize',18)

%%

figure; hold on;
errorbar(T,nanmean(squeeze(mean_R_data(:,:,4)),2),nanstd(squeeze(mean_R_data(:,:,4))'),'r*')
errorbar(T,nanmean(squeeze(mean_R_data(:,:,3)),2),nanstd(squeeze(mean_R_data(:,:,3))'),'k*')
plot(T,nanmean(squeeze(mean_R_data(:,:,4)),2),'r')
plot(T,nanmean(squeeze(mean_R_data(:,:,3)),2),'k')
title('Correlation on data')
ylim([-0.2 0.5])
set(gca,'FontSize',18)

%%
figure; hold on; 
errorbar(T,nanmean(squeeze(mean_R(:,:,1)),2),nanstd(squeeze(mean_R(:,:,1))'),'r*')
errorbar(T,nanmean(squeeze(mean_R(:,:,2)),2),nanstd(squeeze(mean_R(:,:,2))'),'b*')
errorbar(T,nanmean(squeeze(mean_R(:,:,3)),2),nanstd(squeeze(mean_R(:,:,3))'),'k*')
plot(T,nanmean(squeeze(mean_R(:,:,1)),2),'r')
plot(T,nanmean(squeeze(mean_R(:,:,2)),2),'b')
plot(T,nanmean(squeeze(mean_R(:,:,3)),2),'k')
title('Correlation on ising model data')
ylim([-0.2 0.5])

set(gca,'FontSize',18)