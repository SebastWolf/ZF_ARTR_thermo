%% ===========================  Figure 2AC ============================================

% Display histograms of mL, mR from data & ising

clear
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];

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
    
    figure('Name', fn);
    subplot(2,2,1);
    histogram2(mLdata, mRdata, nbins, 'FaceColor', 'flat', ...
        'Normalization', 'pdf');
    ax = gca;
    ax.ZScale = 'log';
    title(['Data ' dirisin(idf).name])
    ax.View = [50, 50];
    xlim([0, 1]);
    ylim([0, 1]);
    grid off
    xlabel('m_L'); ylabel('m_R'); zlabel('pdf');
    ax.XTick = [0, 0.5, 1];
    ax.YTick = [0, 0.5, 1];
    ax.ColorScale = 'log';
    c = colorbar;
    title(c, 'pdf');
    caxis(ax, [1e-3, 10]);
    
    subplot(2,2,2); hold on
    plot(time, mRdata);
    plot(time, mLdata);
    legend('mR', 'mL');
    xlim([400, 1000]);
    ylim([0, 1]);
    xlabel('time (s)');
    ylabel('spikes');
    
    % - ISING
    mLdata = D.Dinference_corr.mleft;
    mRdata = D.Dinference_corr.mright;
    
%     figure;
    subplot(2,2,3);
    histogram2(mLdata, mRdata, nbins, 'FaceColor', 'flat', ...
        'Normalization', 'pdf');
    
    ax = gca;
    ax.ZScale = 'log';
    title(['Ising ' dirisin(idf).name(1:end-4)])
    ax.View = [50, 50];
    xlim([0, 1]);
    ylim([0, 1]);
    grid off
    xlabel('m_L'); ylabel('m_R'); zlabel('pdf');
    ax.XTick = [0, 0.5, 1];
    ax.YTick = [0, 0.5, 1];
    ax.ColorScale = 'log';
    c = colorbar;
    title(c, 'pdf');
    caxis(ax, [1e-3, 10]);
    
    subplot(2,2,4); hold on
    plot(mRdata);
    plot(mLdata);
    legend('mR', 'mL');
    xlim([400, 400+54*600]);
    ylim([0 1]);
    xlabel('time (M.C. round)');
    ylabel('spikes');
    
end