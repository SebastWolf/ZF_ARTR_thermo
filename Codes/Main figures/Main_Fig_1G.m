%% ===========================  Figure 1GH ============================================

% Plot PSD of left/right telegraph signal in behavior

% close all
clear
clc

% --- Parameters
path_princ = '/Users/wolf/Documents/CNRS/Papers/Elife 2023/'
datapath = [path_princ 'Data availability/Thermotaxis/Behavior data/']
filename = [datapath 'allsequences_perbatch.mat']

% * Temperatures
T = [18, 22, 26, 30, 33];

% * Fit parameters
theta_threshold = 10;   % deg

% * PSD parameters
freq = logspace(log10(5e-2), log10(1), 100)';
nw = 4;

% * Figures
colors = rainbow(numel(T));
% --- Prepare figures
figspec = figure; hold on;
figspec.Name = 'spectrumByExp';
fignu = figure; hold on
fignu.Name = 'kflip';
erb = cell(numel(T), 1);

data = load(filename);

% --- Processing
kflip_psd = cell(numel(T), 1);
for idT = 1:numel(T)

    temperature = T(idT);

    % Get corresponding data
    all_turnangle = data.dtheta{idT};
    all_time = data.bouttime{idT};
    framerate = mean([data.framerates{idT}{:}]);
    nexp = numel(all_turnangle);

    % Get PSD from telegraph signal
    kflip_psd{idT} = NaN(nexp, 1);
    apxx = NaN(nexp, numel(freq));
    ftplot = NaN(nexp, numel(freq));

    for idx_exp = 1:nexp

        turnangle = all_turnangle{idx_exp};
        time = all_time{idx_exp};

        nseq = size(turnangle, 2);

        pxx = NaN(nseq, numel(freq));

        for seq = 1:nseq
            seq
            % Get sequence
            dtheta = turnangle(:, seq);
            dtheta(isnan(dtheta)) = [];

            % Trinarize
            dtheta = dtheta - mean(dtheta);  % bias correction
            dthetacopy = dtheta;
            dtheta(abs(dthetacopy) <= theta_threshold) = 0;
            dtheta(dthetacopy > theta_threshold) = 1;
            dtheta(dthetacopy < -theta_threshold) = -1;

            % Reconstruct telegraph time-continuous signal
            dtheta(dtheta == 0) = NaN;      % remove forwards
            t = time(:, seq);
            t(isnan(t)) = [];
            tvec = t(1:end - 1) - t(1);     % start at 0
            fvec = round(tvec*framerate);   % convert in frames
            fr = fvec(1):fvec(end);         % create full time vector
            dthetatime = NaN(size(fr));
            dthetatime(fvec + 1) = dtheta;  % fill with actual values
            dthetatime = fillmissing(dthetatime, 'previous');   % fill missing with previous
            dthetatime = fillmissing(dthetatime, 'next');

            % Get PSD
            pxx(seq, :) = pmtm(dthetatime, nw, freq, framerate);
        end

        % Average PSDs
        mpxx = nanmean(pxx, 1);

        % --- Fit
        % - PXX
        lorenfit = @(nu, I, x) (I*4*nu)./(4*nu^2 + (2*pi*x).^2);
        fo = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0 0], 'Upper', [10 100], 'StartPoint', [0.25 1]);
        mdl = fittype(lorenfit, 'options', fo);
        ft_psd = fit(freq, mpxx', mdl);

        kflip_psd{idT}(idx_exp) = ft_psd.nu;

        apxx(idx_exp, :) = mpxx;
        ftplot(idx_exp, :) = ft_psd(freq)';

    end

    mapxx = mean(apxx, 1);
    eapxx = std(apxx, [], 1)./sqrt(size(apxx, 1));
    mftpl = mean(ftplot, 1);

    errorlogged = (eapxx./mapxx).*10.*log10(mapxx);

    % --- Display
    lineopt = struct;
    lineopt.DisplayName = ['T = ' num2str(temperature) '°C'];
    lineopt.Color = colors(idT, :);
    lineopt.MarkerSize = 6;
    lineopt.Marker = '.';
    lineopt.LineStyle = 'none';
    shadopt = struct;
    shadopt.FaceAlpha = 0.275;

    set(0, 'CurrentFigure', figspec)
    erb{idT} = errorshaded(freq, 10*log10(mapxx), errorlogged, 'line', lineopt, 'patch', shadopt);
    p = plot(freq, 10*log10(mftpl), 'Color', colors(idT, :));
    p.LineWidth = 1;

end

ax = gca; ax.XScale = 'log';
xlabel(ax, 'Frequency [Hz]');
ylabel(ax, '<Power>_{T} [dB/Hz]');
axis(ax, 'square');
legend(ax, [erb{:}], 'Location', 'southwest');
ax.XTick = [0.1, 1];
ylim(ax, [-25, 5]);
pp = plot([NaN NaN], [NaN, NaN], 'k');
pp.DisplayName = 'Fit';

% --- Display
maxsize = max(cellfun(@numel, kflip_psd));
kflip_psd = cellfun(@(x) cat(1, x, NaN(maxsize - numel(x), 1)), kflip_psd, 'UniformOutput', false);
g = cell(size(kflip_psd));
for idx = 1:size(kflip_psd, 1)
    g{idx} = idx.*ones(numel(kflip_psd{idx}), 1);
end
g = cat(1, g{:});
labels = arrayfun(@(x) [num2str(x) '°C'], T, 'UniformOutput', false);
X = cat(1, kflip_psd{:});

set(0, 'CurrentFigure', fignu);
pointopt = struct;
pointopt.Color = colors;
pointopt.jitter = 0;
boxopt = struct;
boxopt.BoxColors = colors;
boxopt.Labels = labels;

beautifulbox(X, g, pointopt, boxopt);

t = title('k_{flip} from PSD fit');
yl = ylabel('k_{flip} (s^{-1})');
axis([-Inf Inf 0 1.1]);
axis square

% --- Save file
save(outfname, 'kflip_psd');