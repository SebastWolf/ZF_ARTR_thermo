%% ===========================  Figure 3A ============================================

% Compare Data & Ising (mL, mR) histograms

clear all
clc
close all

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
hpath_star = @(fish_num) [isinpath 'T*_Fish ' num2str(fish_num) '.mat'];

fish_num = [2 3 4 5 6 7 11 12 13 14 15 16 17];
T = [18, 22, 26, 30, 33];
cmp = [0 0 1; 0 0.5 1; 0 1 0; 1 0.6 0; 1 0 0];
pcount = 0;

% --- Parameters
nbins = 10;
ebins = linspace(0, 1, nbins + 1);
pseudocount = 1;

T = [18, 22, 26, 30, 33];
colors = lines(numel(T));

% --- Init.
dirisin = dir([isinpath '*.mat']);
K = NaN(numel(dirisin), 1);
Krd = NaN(numel(dirisin), numel(dirisin));

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
        
        load([isinpath '/T' num2str(T(idT)) '_Fish ' num2str(fish_num(id_fish)) '.mat'])

        id_fish
        t_files
        
        % load data
       
        list_out_25 = Dinference_corr.out_25;
        data_BM_indpdt = Dinference_corr.synth_indpdt;
        mLdata = Dinference_corr.mleft_data(list_out_25);
        mRdata = Dinference_corr.mright_data(list_out_25);
        mLisin = Dinference_corr.mleft;
        mRisin = Dinference_corr.mright;
        mLindpdt = mean(double(data_BM_indpdt(:,Dinference_corr.L_reg)),2);
        mRindpdt = mean(double(data_BM_indpdt(:,Dinference_corr.R_reg)),2);

        
        % Get counts in each bins
        countdata = histcounts2(mLdata, mRdata, ebins, ebins);
        countisin = histcounts2(mLisin, mRisin, ebins, ebins);
        countindpdt = histcounts2(mLindpdt, mRindpdt, ebins, ebins);

        % Add pseudocounts
        countdata = countdata + pseudocount;
        countisin = countisin + pseudocount;
        countindpdt = countindpdt + pseudocount;

        % Normalize to get pdf
        pddata = countdata./sum(countdata, 'all'); % proba
        pdisin = countisin./sum(countisin, 'all');
        pdindpdt = countindpdt./sum(countindpdt, 'all');
        pddata = pddata(:);
        pdisin = pdisin(:);
        pdindpdt = pdindpdt(:);

        
        
        % Compute KL divergence
        pcount = pcount + 1;
        K(pcount) = sum( pddata.*log10(pddata./pdisin) );
        Kindpdt(pcount) = sum( pddata.*log10(pddata./pdindpdt) );

        % - Data i vs Ising j
        
        for id_fish_rd = 1:numel(fish_num)
            
            dpath_rd = hpath_star(fish_num(id_fish_rd));
            files_rd = dir(dpath_rd);
            
            clear files_T_rd
            pcount_rd = 0;
            for t_files_rd = 1:length(files_rd)
                for idT_rd = 1:numel(T)
                    if isempty(strfind(files_rd(t_files_rd).name,num2str(T(idT_rd))))==0;
                        files_T_rd(t_files_rd) = T(idT_rd);
                    end
                end
                idT_rd = find(files_T_rd(t_files_rd)==T);
                
                Drd = load([isinpath '/T' num2str(T(idT_rd)) '_Fish ' num2str(fish_num(id_fish_rd)) '.mat']);
                
                mLrd = Drd.Dinference_corr.mleft;
                mRrd = Drd.Dinference_corr.mright;
                
                % Get counts
                countrd = histcounts2(mLrd, mRrd, ebins, ebins);
                
                % Add pseudocounts
                countrd = countrd + pseudocount;
                
                % Normalize to get proba
                pdrd = countrd./sum(countrd, 'all');               % proba
                
                pdrd = pdrd(:);
                
                % KL divergence
                pcount_rd = pcount_rd + 1;
                Krd(pcount,pcount_rd) = sum( pddata.*log10(pddata./pdrd) );
                
            end
        end
    end
end

%%

bin_size = 0.05;
vect_edges = [round(min([Kindpdt,Krd(:)',K'])/bin_size)*bin_size:bin_size:round(max([Kindpdt,Krd(:)',K'])/bin_size)*bin_size]

figure('Name', 'KLdiv'); hold on
hd = histogram(K, [vect_edges],'Normalization', 'pdf', ...
    'EdgeColor', 'none', 'FaceColor', [63, 143, 41]./255, 'FaceAlpha', 0.8);
hi = histogram(Krd,  [vect_edges],'Normalization', 'pdf', ...
    'EdgeColor', 'none', 'FaceColor', [222, 25, 35]./255, 'FaceAlpha', 0.8);
hi = histogram(Kindpdt,  [vect_edges],'Normalization', 'pdf', ...
    'EdgeColor', 'none', 'FaceColor', [36, 25, 155]./255, 'FaceAlpha', 0.8);
xlabel('KL divergence');
ylabel('pdf');
legend('Test data_i/Ising_i', 'Test data_i/Ising_{j\neq i}', 'Test data_i/Indpdt_i','Location', 'northeast');
%xlim([0, 2]);
grid off
set(gca,'FontSize',18)
