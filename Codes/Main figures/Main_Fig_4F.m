%% ===========================  Figure 4F ============================================

% Persistence time of the mean-field ARTR model for all ?fish and runs at dfferent experimental 
%temperatures. Each dot refers to one ?Fish at one temperature, colors encode temperature.

clear all
clc
close all
    
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
load([pwd '/Data availability/Mean_field_Parameters.mat'])
hpath_star = @(fish_num) [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J/T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];

% --- Parameters
thresh = 0.1;
thresh_MF = 0.1;
prctage = 50;       % percentage of activity as threshold
nrmd = 1;           % remove persistence < nrmd seconds in data
nrmi = 30;         % remove persistence < nrmi steps in ising
nbins = 50;         % bins for persistence time distributions
ebins = linspace(0, 20, nbins);

dispall = 'n';
colors = lines(2);

tcolors = lines(numel(T));


pcount = 1;

for id_fish = 1:numel(fish_num)
    
    dpath = hpath_star(fish_num(id_fish));
    clear files_T
    files = dir(dpath);
    
    for t_files = 1:length(files)
        for idT = 1:numel(T)
            if isempty(strfind(files(t_files).name,num2str(T(idT))))==0;
                files_T(t_files) = T(idT);
            end
        end
        idT = find(files_T(t_files)==T);
        
        line_selected = find(table_mean_field(:,1)==T(idT) & table_mean_field(:,2)==fish_num(id_fish));
        
        load([files(t_files).folder filesep files(t_files).name]);
               
        temperatures(pcount) = T(idT);

        L_reg = Dinference_corr.L_reg;
        R_reg = Dinference_corr.R_reg;
        N_L = length(L_reg);
        N_R = length(R_reg);        
        J_inter = Dinference_corr.J(1:N_L,N_L+1:N_L+N_R);
        J_left = Dinference_corr.J(1:N_L,1:N_L);
        J_right = Dinference_corr.J(N_L+1:N_L+N_R,N_L+1:N_L+N_R) ;
               
        M_h_R = table_mean_field(line_selected,7);
        M_h_L = table_mean_field(line_selected,6);;
        M_J_R = nanmean(J_right(:));
        M_J_L = nanmean(J_left(:));        
        M_I_I = nanmean(J_inter(:));
        
        J_J_L = table_mean_field(line_selected,3);
        J_J_R = table_mean_field(line_selected,4);
        I_I = table_mean_field(line_selected,5);
        
        % activity
                        
        N_L_noise = table_mean_field(line_selected,8);
        N_R_noise = table_mean_field(line_selected,9);
        N_L = N_L_noise;
        N_R = N_R_noise;
        
        
        
        tau = 1;
            dt = 0.001;
            
            thr_lim = 1e-4;
            
            sigma_L = sqrt(tau/dt)*sqrt(2)
            sigma_R = sqrt(tau/dt)*sqrt(2);
            T_time = 10000;
            mr = thr_lim*ones(1,T_time);
            ml = thr_lim*ones(1,T_time);
            
            N_L = N_L_noise;
            N_R = N_R_noise;
            
            
            for i = 1:T_time-1
                
                %ml
                xi_L = sigma_L*randn(1);
                xi_L_note(i) = xi_L;
                ml(i+1) = ml(i)-dt*(1/tau)*(RFIM_free_energy_DM_U_paper_2(ml(i),mr(i),N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I)+xi_L);
                if ml(i+1)<=thr_lim; ml(i+1) = thr_lim; end
                if ml(i+1)>=1-thr_lim; ml(i+1) = 1-thr_lim; end
                %RFIM_free_energy_DM_U_paper(ml(i),mr(i),N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I)
                %mr
                xi_R = sigma_R*randn(1);
                mr(i+1) = mr(i)-dt*(1/tau)*(RFIM_free_energy_DM_V_paper_2(ml(i),mr(i),N_L,N_R,M_h_L,M_h_R,M_J_L,M_J_R,M_I_I,J_J_L,J_J_R,I_I)+xi_R);
                if mr(i+1)<=thr_lim; mr(i+1) = thr_lim; end
                if mr(i+1)>=1-thr_lim; mr(i+1) = 1-thr_lim; end
                
            end
    
%%       
% - DATA
    mLdata = Dinference_corr.mleft_data;
    mRdata = Dinference_corr.mright_data;
    time = Dinference_corr.time;
    
    framerate = 1./mean(diff(time));
    time = linspace(time(1), time(end), numel(mLdata));
    nrm = nrmd*framerate;
    
    % Get threshold
    threshLd = thresh;
    threshRd = thresh;
%     threshLd = prctile(mLdata, prctage);
%     threshRd = prctile(mRdata, prctage);
    
%     ths(idf, :) = [threshLd, threshRd];
    
    % Get logical
    logiLd = mLdata > threshLd;
    logiRd = mRdata > threshRd;
    
    % Get persistence time
    LL = bwconncomp(logiLd);
    LR = bwconncomp(logiRd);
    
    persisLd = cellfun(@numel, LL.PixelIdxList);
    persisRd = cellfun(@numel, LR.PixelIdxList);
    
    % Find too short periods
    tormLd = persisLd <= nrm;
    tormRd = persisRd <= nrm;
    
    % Filter them out
    tormindLd = LL.PixelIdxList;
    tormindLd(~tormLd) = [];
    tormindLd = cat(1, tormindLd{:});
    tormindRd = LR.PixelIdxList;
    tormindRd(~tormRd) = [];
    tormindRd = cat(1, tormindRd{:});
    logiLd(tormindLd) = 0;
    logiRd(tormindRd) = 0;
    
    persisLd(tormLd) = [];
    persisRd(tormRd) = [];
    
    % Convert to time
    persisLd = persisLd'./framerate;
    persisRd = persisRd'./framerate;
    
    % Store
    pool = cat(1, persisLd, persisRd);
    dpperstime{pcount} = pool;
    %dplabels{pcount} = repmat(fn, numel(pool), 1);
    dptemperature{pcount} = repmat(T, numel(pool), 1);
    mdata(pcount) = mean(pool);            
    stable_data_all_fit(idT,id_fish) =  mean(pool);

            
%%            
   % - MEAN FIELD
    mLisin = ml;
    mRisin = mr;
    
    % Get threshold
    threshLi = thresh_MF;
    threshRi = thresh_MF;
    
    % Get logical
    logiLi = mLisin > threshLi;
    logiRi = mRisin > threshRi;
    
    % Get persistence time
    LL = bwconncomp(logiLi);
    LR = bwconncomp(logiRi);
    
    persisLi = cellfun(@numel, LL.PixelIdxList);
    persisRi = cellfun(@numel, LR.PixelIdxList);
    
    % Find too short periods
    tormLi = persisLi <= nrmi;
    tormRi = persisRi <= nrmi;
    
    % Filter them out
    tormindLi = LL.PixelIdxList;
    tormindLi(~tormLi) = [];
    tormindLi = cat(1, tormindLi{:});
    tormindRi = LR.PixelIdxList;
    tormindRi(~tormRi) = [];
    tormindRi = cat(1, tormindRi{:});
    logiLi(tormindLi) = 0;
    logiRi(tormindRi) = 0;
    
    persisLi(tormLi) = [];
    persisRi(tormRi) = [];
    
    % Convert to rounds
    persisLi = persisLi'./framerate;
    persisRi = persisRi'./framerate;
    
    % Store
    pool = cat(1, persisLi, persisRi);
    ipperstime{pcount} = pool;
    %iplabels{pcount} = repmat(fn, numel(pool), 1);
    iptemperature{pcount} = repmat(T, numel(pool), 1);
    misin(pcount) = mean(pool);    
    stable_all_fit(idT,id_fish) =  mean(pool);
    
        pcount = pcount +1
        
    end
end
%%
figure; 

cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

for idT = 1:numel(T)
    
figure(100); hold on
    plot(stable_data_all_fit(idT,:),stable_all_fit(idT,:),'.','MarkerSize',20, 'MarkerEdgeColor',cmp(idT,:))
end


%xlim([0 30])
%ylim([0 100])

%set(gca,'YScale','log')
set(gca,'FontSize',18)
xlabel('Data')
ylabel('Model')


stable_all_fit(stable_all_fit==0) = NaN;
stable_data_all_fit(stable_data_all_fit==0) = NaN;

xx = stable_all_fit(:);
yy = stable_data_all_fit(:);

xx(isnan(xx)==1) = [];
yy(isnan(yy)==1) = [];
figure; plot(stable_data_all_fit,stable_all_fit,'k.')
corrcoef(xx,yy)


