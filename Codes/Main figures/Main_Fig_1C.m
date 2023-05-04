%% ===========================  Figure 1C ============================================
% DSP on ARTR, pool per temperature

% close all
clear
clc

% --- Definitions
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
datapath = [pwd 'Data availability/Thermotaxis/Behavior data/ARTR_export/']

T = [18, 22, 26, 30, 33];

% --- Parameters

% Use the same time vector for all experiment
lag = 10;               % discard fist lag seconds
fs = 5;                 % reference sample rate (Hz)
tstop = 1100;           % Last time to take (s)
reftime_base = 0:1/fs:tstop-1/fs;

% PSD
startf = 0.5*10^-2;
stopf = 0.8;
npoints = 150;          % number of frequency points
freq = logspace(log10(startf), log10(stopf), npoints)';
nw = 4;

% Figure
fig = figure; hold on
colors = rainbow(5);
plt0 = cell(numel(T), 1);
offset = linspace(0, 40, numel(T));    % offset in dB for each curve

% --- Main
nu = NaN(numel(T), 1);
conf_int = NaN(numel(T), 2);
for idT = 1:numel(T)
    
    fprintf('Computing PSD for T = %i°C', T(idT)); tic
    
    % Get experiment list
    list = dir([datapath 'T' num2str(T(idT)) '_Fish*.mat']);

    pxx = NaN(numel(freq), numel(list));
    for idrun = 1:numel(list)
        
        data = load([datapath list(idrun).name]);
        
        mright = data.mright;
        mleft = data.mleft;
        time = data.time;
        motions = data.motions;
        framerate = 1/mean(diff(time));
        
        mright(motions) = NaN;
        mleft(motions) = NaN;
        fstart = fix(lag*framerate);
        fstop = fix(tstop*framerate);
        time = time(fstart:fstop);
        mright = mright(fstart:fstop);
        mleft = mleft(fstart:fstop);
        
        % - Rescale
        mright = rescale(mright, 0, 1);
        mleft = rescale(mleft, 0, 1);
        mright(isnan(mright)) = 0;
        mleft(isnan(mleft)) = 0;
        
        % - Diff
        artr = mright - mleft;
        time = double(time);
        time = time - time(1);
        artr = interp1(time, artr, reftime_base, 'pchip', NaN);
        reftime = reftime_base;
        reftime(isnan(artr)) = [];
        artr(isnan(artr)) = [];
        artr = rescale(artr, -1, 1);
        artr = double(artr - mean(artr));
        
        artr = artr';
        time = time';
   
        artr(artr < 0) = 0;
        artr(artr > 0) = 1;
        
        % - PSD
        pxx(:, idrun) = pmtm(artr, [], freq, fs);
       
        fprintf('.');
        
    end
    
    % - Mean PSD
    mpxx = nanmean(pxx, 2);
    epxx = nanstd(pxx, [], 2)./sqrt(size(pxx, 2));
    
    % - Fit
    fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0], 'Upper', [20 Inf], 'StartPoint', [0.5 0.5]);
    lorenfit = @(nu, I, x) 10*log10((I*4*nu)./(4*nu^2 + (2*pi*x).^2));
    mdl = fittype(lorenfit, 'options', fo);
    ft = fit(freq, 10*log10(mpxx), mdl);
     
    nu(idT) = ft.nu;
    dummy = confint(ft);
    conf_int(idT, :) = dummy(:, 1)';
    
    fprintf('\tDone (%2.2fs).\n', toc);
    
    % --- Display
    set(0, 'CurrentFigure', fig);
    lineopt = struct;
    lineopt.DisplayName = ['T = ' num2str(T(idT)) '°C'];
    lineopt.Color = colors(idT, :);
    lineopt.MarkerSize = 6;
    lineopt.LineStyle = 'none';
    shadopt = struct;
    shadopt.FaceAlpha = 0.275;
    plterror = (epxx./mpxx).*10.*log10(mpxx);
    pltmpxx = 10*log10(mpxx) + offset(idT);
    
    plt0{idT} = errorshaded(freq, pltmpxx, plterror, 'line', lineopt, 'patch', shadopt);
    
    plt1 = plot(freq, ft(freq) + offset(idT));

    plt1.Color = colors(idT, :);
    
    plt2 = plot([ft.nu/(2*pi), ft.nu/(2*pi)], [ft(ft.nu/(2*pi))-5+offset(idT), ft(ft.nu/(2*pi))+5+offset(idT)]);
    plt2.Color = colors(idT, :);
    plt2.LineStyle = '--';
    plt2.LineWidth = 1;
end

% - Cosmetic
ax = gca;
ax.XScale = 'log';
xlabel('frequency (Hz)');
ylabel('<power>_T (dB/Hz)');
fkplt = plot([NaN NaN], [NaN NaN], 'k');
fkplt.DisplayName = 'Fit';
leg = legend([plt0{:}, fkplt]);
leg.Location = 'southwest';
axis square
ylim(ax, [-35, 45]);
fig.Name = 'PooledPSD';