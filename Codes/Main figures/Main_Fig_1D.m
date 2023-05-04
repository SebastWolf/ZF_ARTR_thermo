%% ===========================  Figure 1D ============================================
% DSP on ARTR neurons, extract nu parameter for each fish

% close all
clear
clc

% --- Definitions
path_princ = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
datapath = [path_princ 'Data availability/Thermotaxis/Behavior data/ARTR_export/']
outfname = [path_princ 'Data availability/artr_nu_psd.mat'];

T = [18, 22, 26, 30, 33];

% --- Parameters

% Use the same time vector for all experiment
lag = 10;               % discard fist lag seconds
fs = 5;                 % reference sample rate (Hz)
tstop = 1100;           % Last time to take (s)
reftime = 0:1/fs:tstop-1/fs;

% PSD
startf = 0.5*10^-2;
stopf = 0.8;
npoints = 150;          % number of frequency points
freq = logspace(log10(startf), log10(stopf), npoints)';
nw = 4;

% Figure
fignu = figure('Name', 'NuFromPSD');
axnu = axes(fignu); hold(axnu, 'on');
colors = rainbow(5);

% Preparation
nu = cell(numel(T), 1);
fid = cell(numel(T), 1);
nufish = cell(50, 1);

% --- Main
for idT = 1:numel(T)
    
    fprintf('Computing PSD for T = %i°C', T(idT)); tic
    
    % Get experiment list
    list = dir([datapath 'T' num2str(T(idT)) '_Fish*.mat']);
    
    nu{idT} = NaN(numel(list), 1);
    for idrun = 1:numel(list)
                
        data = load([datapath list(idrun).name]);
        fishid = str2double(list(idrun).name(9:10));
        
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
        hbo = mright - mleft;
        time = double(time);
        time = time - time(1);
        hbo = interp1(time, hbo, reftime, 'pchip', NaN);
        hbo(isnan(hbo)) = [];
        hbo = rescale(hbo, -1, 1);
        hbo = double(hbo - mean(hbo));
        
        hbo = hbo';

        hbo(hbo < 0) = 0;
        hbo(hbo > 0) = 1;

        % - PSD
        pxx = pmtm(hbo, nw, freq, fs);

        % - Fit
        fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0], 'Upper', [Inf Inf], 'StartPoint', [1 1/30]);
        
        psdfit_o = @(I, nu, x) (I*4*nu)./(4*nu^2 + (2*pi*x).^2);
        mdl_o = fittype(@(I, nu, x) 10*log10(psdfit_o(I, nu, x)), 'options', fo);
        ft_o = fit(freq, 10*log10(pxx), mdl_o);

        % - Store
        nu{idT}(idrun) = ft_o.nu;
        nufish{fishid}(idT) = ft_o.nu;
        
        fprintf('.');
        
    end
    
    fprintf('\t Done (%2.2fs).\n', toc);
end

% --- Display
pointopt = struct;
pointopt.Color = colors;
pointopt.jitter = 0;
boxopt = struct;
boxopt.BoxColors = colors;
labels = arrayfun(@(x) [num2str(x) '°C'], T, 'UniformOutput', false);
boxopt.Labels = labels;
boxopt.Notch = 'off';

% * nu
maxsize = max(cellfun(@numel, nu));
nu = cellfun(@(x) cat(1, x, NaN(maxsize - numel(x), 1)), nu, 'UniformOutput', false);
g = cell(size(nu));
for idx = 1:size(nu, 1)
    g{idx} = idx.*ones(numel(nu{idx}), 1);
end
g = cat(1, g{:});

X = cat(1, nu{:});

set(0, 'CurrentFigure', fignu);
set(fignu, 'CurrentAxes', axnu);
hold(axnu, 'on');

beautifulbox(X, g, pointopt, boxopt);
title('Rate parameter');
ylabel('\nu (s^{-1})');
axis([-Inf Inf 0 1.1]);
axis(axnu, 'square');
ax = gca;
ax.YGrid = 'off';

% plot same fish lines
for id = 1:size(nufish, 1)
    
    fnu = nufish{id};
    
    x = find(fnu);
    
    p = plot(axnu, x, fnu(x));
    
    if ~isempty(p)
        p.Color = [0, 0, 0, 0.25];
        p.LineWidth = 1;
    end
end

save(outfname, 'nu');