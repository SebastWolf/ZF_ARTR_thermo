%% ================================ Figure 5E ============================================
% An example Monte Carlo activity trace generated with the modified Ising model.

clear all
close all
clc

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
hpath_light = [pwd filesep 'Data availability' filesep 'Light stim' filesep ];

dpath = hpath_light;
files = dir([dpath '/*.mat']);
nfiles = numel(files);

list_idf = [1:numel(files)]
idf = list_idf(3)

load([files(idf).folder filesep files(idf).name])
CONFIG_h_sup = Dinference.CONFIG_under_stim_LONG;
s_plus = Dinference.s_plus;
s_moins = Dinference.s_moins;

L_vec = Dinference.L_reg;
R_vec = Dinference.R_reg;
mL_MC = mean(CONFIG_h_sup(:,L_vec),2);
mR_MC = mean(CONFIG_h_sup(:,R_vec),2);

vect_stim_L = 0*mL_MC;
vect_stim_R = 0*mL_MC;
vect_stim_L(s_plus) = 1;
vect_stim_R(s_moins) = 1;

figure; hold on
plot(mL_MC,'r','LineWidth',2)
plot(mR_MC,'b','LineWidth',2)
plot(vect_stim_L,'r')
plot(vect_stim_R,'b')
xlabel('MC');
ylabel('m_{L,R}')
set(gca,'FontSize',18)

xlim([0 1000])
