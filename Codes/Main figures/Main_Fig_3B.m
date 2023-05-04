%% ===========================  Figure 3B ============================================

% Get persistence time from each side of ARTR

clear all
close all
clc

% --- Definitions
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
fnout = [pwd filesep 'Data' filesep 'Thermotaxis' filesep 'Matfiles' filesep 'persistime.mat'];

allT = [18, 22, 26, 30, 33];

% --- Parameters
thresh = 0.1;
prctage = 50;       % percentage of activity as threshold
nrmd = 1;           % remove persistence < nrmd seconds in data
nrmi = 30;         % remove persistence < nrmi steps in ising
nbins = 50;         % bins for persistence time distributions
ebins = linspace(0, 20, nbins);

dispall = 'n';
colors = lines(2);

%tcolors = tempcol(numel(allT));
tcolors = lines(numel(allT));

% --- Init.
dirisin = dir([isinpath '*.mat']);

dpperstime = cell(numel(dirisin), 1);
dplabels = cell(numel(dirisin), 1);
dptemperature = cell(numel(dirisin), 1);
ipperstime = cell(numel(dirisin), 1);
iplabels = cell(numel(dirisin), 1);
iptemperature = cell(numel(dirisin), 1);
fishids = NaN(numel(dirisin), 1);
temperatures = NaN(numel(dirisin), 1);
mdata = NaN(numel(dirisin), 1);
misin = NaN(numel(dirisin), 1);
ths = NaN(numel(dirisin), 2);



for idf = 1:numel(dirisin)
    idf
    % Get data
    fnisin = [isinpath dirisin(idf).name];
    
    % Load data
    D = load(fnisin);
    
    % Meta
    loc = strfind(dirisin(idf).name, '.mat');
    idata = str2double(dirisin(idf).name(loc-2:loc-1));
    fn = [dirisin(idf).name(1:3) '-Fish' num2str(idata, '%02i')];
    T = D.Dinference_corr.T;

    temperatures(idf) = T;
    fishids(idf) = idata;
    
    % - DATA
    time = D.Dinference_corr.time;
    mLdata = D.Dinference_corr.mleft_data;
    mRdata = D.Dinference_corr.mright_data;
    
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
    dpperstime{idf} = pool;
    dplabels{idf} = repmat(fn, numel(pool), 1);
    dptemperature{idf} = repmat(T, numel(pool), 1);
    mdata(idf) = mean(pool);
    
    % - ISING
    mLisin = D.Dinference_corr.mleft;
    mRisin = D.Dinference_corr.mright;
    
    % Get threshold
    threshLi = thresh;
    threshRi = thresh;
%     threshLi = prctile(mLisin, prctage);
%     threshRi = prctile(mRisin, prctage);
    
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
    ipperstime{idf} = pool;
    iplabels{idf} = repmat(fn, numel(pool), 1);
    iptemperature{idf} = repmat(T, numel(pool), 1);
    misin(idf) = mean(pool);
    
    % - DISPLAY
    switch dispall
        case 'y'
            figure('Name', fn);
            
            subplot(2, 2, 1); hold on; title('Data');
            A = area(time, logiRd);
            A.EdgeColor = 'none';
            A.FaceColor = colors(1, :);
            A.FaceAlpha = 0.4;
            plot([time(1), time(end)], [threshRd, threshRd], 'Color', colors(1, :));
            tmp = mRdata;
            tmp(logiRd) = NaN;
            p = plot(time, tmp, 'Color', colors(1, :));
            p.Color(4) = 0.2;
            tmp = mRdata;
            tmp(~logiRd) = NaN;
            plot(time, tmp, 'Color', colors(1, :));
            axis tight
            grid off
            xlabel('time (s)');
            ylabel('m_R (spikes)');
            
            subplot(2, 2, 3); hold on
            A = area(time, logiLd);
            A.EdgeColor = 'none';
            A.FaceColor = colors(2, :);
            A.FaceAlpha = 0.4;
            plot([time(1), time(end)], [threshLd, threshLd], 'Color', colors(2, :));
            tmp = mLdata;
            tmp(logiLd) = NaN;
            p = plot(time, tmp, 'Color', colors(2, :));
            p.Color(4) = 0.25;
            tmp = mLdata;
            tmp(~logiLd) = NaN;
            plot(time, tmp, 'Color', colors(2, :));
            axis tight
            grid off
            xlabel('time (s)');
            ylabel('m_L (spikes)');
            
            subplot(2, 2, 2); hold on; title('Ising')
            A = area(1:numel(mRisin), logiRi);
            A.EdgeColor = 'none';
            A.FaceColor = colors(1, :);
            A.FaceAlpha = 0.4;
            plot([1, numel(mRisin)], [threshRi, threshRi], 'Color', colors(1, :));
            tmp = mRisin;
            tmp(logiRi) = NaN;
            p = plot(tmp, 'Color', colors(1, :));
            p.Color(4) = 0.25;
            tmp = mRisin;
            tmp(~logiRi) = NaN;
            plot(tmp, 'Color', colors(1, :));
            axis tight
            grid off
            xlabel('time (?)');
            ylabel('m_R (spikes)');
            
            subplot(2, 2, 4); hold on;
            A = area(1:numel(mLisin), logiLi);
            A.EdgeColor = 'none';
            A.FaceColor = colors(2, :);
            A.FaceAlpha = 0.4;
            plot([1, numel(mLisin)], [threshLi, threshLi], 'Color', colors(2, :));
            tmp = mLisin;
            tmp(logiLi) = NaN;
            p = plot(tmp, 'Color', colors(2, :));
            p.Color(4) = 0.25;
            tmp = mLisin;
            tmp(~logiLi) = NaN;
            plot(tmp, 'Color', colors(2, :));
            axis tight
            grid off
            xlabel('time (?)');
            ylabel('m_L (spikes)');
            
            sgtitle(fn);
            
    end
    
end
%%
% Pool
dptemperature = cat(1, dptemperature{:});
dplabels = categorical(cellstr(cat(1, dplabels{:})));
dpperstime = cat(1, dpperstime{:});
iptemperature = cat(1, iptemperature{:});
iplabels = categorical(cellstr(cat(1, iplabels{:})));
ipperstime = cat(1, ipperstime{:});
unids = unique(fishids);

temperatures(find(temperatures==27)) = 26;
%%
figure;
[~, I] = ismember(temperatures, allT);
cols = tcolors(I, :);

% - Data
subplot(1, 2, 1); hold on
scatter(temperatures, mdata, [], cols, 'filled');
xlabel('temperature (°C)');
ylabel('<persistence time> (s)');
title('Data');

% add lines
for id = 1:numel(unids)
    
    X = temperatures(fishids == unids(id));
    Y = mdata(fishids == unids(id));
    
    plot(X, Y, 'Color', [0.5, 0.5, 0.5, 0.75]);
    
end

ax = gca;
ax.XTick = allT;
ax.YGrid = 'off';
xlim([17, 34]);

% - Ising
subplot(1, 2, 2); hold on
scatter(temperatures, misin, [], cols, 'filled');
xlabel('temperature (°C)');
ylabel('<persistence time> (M.C. round)');
title('Ising');

% add lines
for id = 1:numel(unids)
    
    X = temperatures(fishids == unids(id));
    Y = misin(fishids == unids(id));
    
    plot(X, Y, 'Color', [0.5, 0.5, 0.5, 0.75]);
    
end

ax = gca;
ax.XTick = allT;
ax.YGrid = 'off';
xlim([17, 34]);

%%
figure('Name', 'Ising_data_persistence'); hold on

% * Conversion time second <-> MC round
ft = fit(mdata, misin, 'a*x');
plot(0:20, ft(0:20), 'k--');

% * Data vs Ising

scatter(mdata, misin, [], cols, 'filled');
xlabel('<persistence time>_{data} (s)');
ylabel('<persistence timen>_{ising} (M.C. rounds)');
xlim([0, 20]);
ylim([0, 180]);
R = corrcoef([mdata, misin]);
text(10, 40, ['R_{Pearson} = ' num2str(R(1, 2), '%0.2f')], ...
    'FontSize', 16);
axis square
title(['threshold = ' num2str(thresh)]);
grid off
