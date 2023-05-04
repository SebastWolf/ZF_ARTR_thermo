%% ================================ SI Figure 3EF ============================================

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];

CCC_ALL = [];
CCC_model_ALL = [];
CCC_no_diag_ALL = [];
CCC_no_diag_model_ALL = [];
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
        load([isinpath '/T' num2str(T(idT)) '_Fish ' num2str(fish_num(id_fish)) '.mat'])

        A_cor_rh23_75 = Dinference_corr.data_rh23_75;
        A_cor_rh23_25 = Dinference_corr.data_rh23_25;
        list_out_25 = Dinference_corr.out_25;
        data_BM_indpdt = double(Dinference_corr.synth_indpdt);
                
        firing_rates_75 = mean(A_cor_rh23_75,1);
        firing_rates__data_BM_75 = mean(data_BM_indpdt,1);
        Lmax = 50000;
        data_indpdt_model = zeros(Lmax,length(firing_rates_75));
        for t=1:Lmax
            data_indpdt_model(t,:) = (rand(length(firing_rates_75),1)-firing_rates_75'<0)';
        end
        
        
        
        h_field = log(firing_rates_75./(1-firing_rates_75));
        proba = proba_indep_model(data_BM_indpdt', h_field');
        LL_Id = log(proba)/length(firing_rates_75);

        pcount = pcount + 1;
            
            figure(2);
            %clf
            hold on;
            plot(nanmean(A_cor_rh23_25,1),nanmean(data_indpdt_model,1),'k.')
            xlabel('Data')
            ylabel('Independent model')
            set(gca,'FontSize',18)
            
            figure(1);
            hold on;
            plot(cov(A_cor_rh23_25),cov(data_BM_indpdt),'k.')
            xlabel('Data')
            ylabel('Independent model')
            set(gca,'FontSize',18)
            
            CCC = cov(A_cor_rh23_25);
            CCC_ALL = [CCC_ALL CCC(:)'];
            
            CCC_model = cov(data_BM_indpdt);
            CCC_model_ALL = [CCC_model_ALL CCC_model(:)'];

            CCC = cov(A_cor_rh23_25);
            CCC_model = cov(data_BM_indpdt);
            for i = 1:length(firing_rates_75); CCC(i,i) = NaN; CCC_model(i,i) = NaN; end
            CCC_no_diag_ALL = [CCC_no_diag_ALL CCC(:)'];            
            CCC_no_diag_model_ALL = [CCC_no_diag_model_ALL CCC_model(:)'];

            figure(3);
            hold on;
            plot(CCC_no_diag_ALL,CCC_no_diag_model_ALL,'k.')
            xlabel('Data')
            ylabel('Independent model')
            set(gca,'FontSize',18)
            
            FFF = nanmean(A_cor_rh23_25,1);
            FFF_ALL = [FFF_ALL FFF];
            
            FFF_model = nanmean(data_BM_indpdt,1);
            FFF_model_ALL = [FFF_model_ALL FFF_model];
            
                
    end
    
end


%%


xmin = -0.3; xmax = 0.4;

fitttCC = fit(CCC_ALL',CCC_model_ALL','poly1')
figure(1);
plot([xmin xmax],fitttCC([xmin xmax]),'r')
xlim([xmin xmax])
ylim([xmin xmax])

%%

xmin = -0.3; xmax = 0.4;
CCC_no_diag_ALL(isnan(CCC_no_diag_ALL)==1) = [];
CCC_no_diag_model_ALL(isnan(CCC_no_diag_model_ALL)==1) = [];

fitttCC = fit(CCC_no_diag_ALL',CCC_no_diag_model_ALL','poly1')
figure(3);
plot([xmin xmax],fitttCC([xmin xmax]),'r')
xlim([xmin xmax])
ylim([-0.05 0.05])

%%

fitttFF = fit(FFF_ALL',FFF_model_ALL','poly1')
xmin = 0; xmax = 1;

figure(2);
plot([xmin xmax],fitttFF([xmin xmax]),'r')

xlim([xmin xmax])
ylim([xmin xmax])