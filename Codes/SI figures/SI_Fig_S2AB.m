% Distribution (1D) of mL and mR

clear all
close all
clc

% --- Definitions

pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];
fnout = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep '/Behavior data/meanmLmR.mat'];


allT = [18, 22, 26, 30, 33];
allids = [2:5, 6, 7, 11:17];

% --- Parameters
nbins = 50;
ebins = linspace(0, 1, nbins + 1);
nbins2D = 10;

dispall = 'n';
colors = lines(2);
tcols = rainbow(numel(allT));

fd = figure('Name', 'DatamDistByT');
axd = gca; hold(axd, 'on');
fi = figure('Name', 'IsingmDistByT');
axi = gca; hold(axi, 'on');
fm = figure('Name', 'MeanDataIsing');
axm = gca; hold(axm, 'on');

% --- Init.
dirisin = dir([isinpath '*.mat']);

am1pdfd = NaN(1, numel(dirisin));
am1pdfi = NaN(1, numel(dirisin));
allpdfd = NaN(nbins, numel(dirisin));
allpdfi = NaN(nbins, numel(dirisin));
fishids = NaN(numel(dirisin), 1);
temperatures = NaN(numel(dirisin), 1);

for idf = 1:numel(dirisin)
    
    fnisin = [isinpath dirisin(idf).name];
    D = load(fnisin);
    
    % Meta
    loc = strfind(dirisin(idf).name, '.mat');
    idata = str2double(dirisin(idf).name(loc-2:loc-1));
    fn = [dirisin(idf).name(1:3) '-Fish' num2str(idata, '%02i')];
    T = D.Dinference_corr.T;
    
    % - DATA
    mLdata = D.Dinference_corr.mleft_data;
    mRdata = D.Dinference_corr.mright_data;
    
    % Get 1D pdf
    [pdfdata, cbins] = computePDF(ebins, [mLdata; mRdata], 'mode', 'edges');
    
    % Get first moment
    am1pdfd(idf) = trapz(cbins, cbins'.*pdfdata);

    % - ISING
    mLisin = D.Dinference_corr.mleft;
    mRisin = D.Dinference_corr.mright;
    
    % Get 1D pdf
    pdfisin = computePDF(ebins, [mRisin; mRdata], 'mode', 'edges');
    
    % Get first moment
    am1pdfi(idf) = mean([mLisin; mRisin]);
    
    % Store
    fishids(idf) = idata;
    temperatures(idf) = T;
    allpdfd(:, idf) = pdfdata;
    allpdfi(:, idf) = pdfisin;
    
    % - DISPLAY
    switch dispall
        case 'y'
            figure('Name', fn);
            subplot(2, 1, 1); hold on
            plot(cbins, pdfdata);
            ax = gca;
            ax.YScale = 'log';
            xlabel('spikes');
            ylabel('pdf');
            title('Data');
            
            subplot(2, 1, 2); hold on
            plot(cbins, pdfisin);
            ax = gca;
            ax.YScale = 'log';
            xlabel('spikes');
            ylabel('pdf');
            title('Ising');
            
            sgtitle(fn);
    end
    
    if T==27; T=26; end
    
    p = plot(axd, cbins, pdfdata);
    p.Color = tcols(ismember(allT, T), :);
    p.Color(4) = 0.2;
    
    p = plot(axi, cbins, pdfisin);
    p.Color = tcols(ismember(allT, T), :);
    p.Color(4) = 0.2;
    
    r = corrcoef([mLisin, mRisin]);
    Rising(idf) = r(2,1);
    r = corrcoef([mLdata, mRdata]);
    Rdata(idf) = r(2,1);
    
end



% Mean over temperature
Tm1pdfd = NaN(numel(allT), 1);
Ts1pdfd = NaN(numel(allT), 1);
Tm1pdfi = NaN(numel(allT), 1);
Ts1pdfi = NaN(numel(allT), 1);

for idT = 1:numel(allT)
    
    % - Mean first moment by temperature
    subm1d = am1pdfd(temperatures == allT(idT));
    Tm1pdfd(idT) = mean(subm1d);
    Ts1pdfd(idT) = std(subm1d)./sqrt(numel(subm1d));
    subm1i = am1pdfi(temperatures == allT(idT));
    Tm1pdfi(idT) = mean(subm1i);
    Ts1pdfi(idT) = std(subm1d)./sqrt(numel(subm1i));
    
    % - Mean pdf by temperature
    subpdfd = allpdfd(:, temperatures == allT(idT));
    subpdfi = allpdfi(:, temperatures == allT(idT));
    
    mpdfd = mean(subpdfd, 2);
    epdfd = std(subpdfd, [], 2)./sqrt(size(subpdfd, 2));
    
    mpdfi = mean(subpdfi, 2);
    epdfi = std(subpdfi, [], 2)./sqrt(size(subpdfi, 2));
        
    plot(axd, cbins, mpdfd, 'Color', tcols(idT, :), 'LineWidth', 2);
    plot(axi, cbins, mpdfi, 'Color', tcols(idT, :), 'LineWidth', 2);
    
    scatter(axm, subm1d, subm1i, 64, tcols(idT, :), 'filled');
end

% Cosmetic
blank = NaN(numel(allT), 2);
dn = arrayfun(@num2str, allT', 'UniformOutput', false);
s = arrayfun(@(x,y,z1,z2,z3) scatter(x, y, [], [z1, z2, z3], 'filled'), ...
    blank(:, 1), blank(:, 2), tcols(:, 1), tcols(:, 2), tcols(:, 3));

grid(axm, 'off');
axis(axm, 'square');
xlim(axm, [0, 0.4]);
ylim(axm, [0, 0.4]);
pp = plot([0, 1], [0, 1], 'k--');
dn{end + 1} = 'y = x';
lg = legend([s(:); pp], dn, 'Location', 'northwest');
title(lg, 'temperature (°C)');
xlabel(axm, 'data <m_{L,R}>');
ylabel(axm, 'Ising <m_{L,R}>');
R = corrcoef([am1pdfd', am1pdfi']);
text(0.25, 0.1, ['R = ' num2str(R(2, 1), '%0.2f')], 'FontSize', 12);
lg.AutoUpdate = 'off';

p = cell(numel(allT), 1);
for idT = 1:numel(allT)
    p{idT} = plot(axi, [NaN, NaN], [NaN, NaN], 'Color', tcols(idT, :));
    p{idT}.DisplayName = num2str(allT(idT));
end

xlabel(axd, 'm_{L&R}');
ylabel(axd, 'pdf');
l = legend(axd, [p{:}]);
title(l, 'temperature (°C)');
axd.YScale = 'log';
title(axd, 'Data');
xlim(axd, [0, 0.3]);
% ylim(axd, [5e-3, 5e-1]);
axis(axd, 'square');

xlabel(axi, 'm_{L&R}');
ylabel(axi, 'pdf');
l = legend(axi, [p{:}]);
title(l, 'temperature (°C)');
axi.YScale = 'log';
title(axi, 'Ising');
xlim(axi, [0, 0.3])

% First moments
figure; hold on
for idT = 1:numel(allT)
    
    % - Mean first moment by temperature
    subm1d = am1pdfd(temperatures == allT(idT));
    subT = repmat(allT(idT), numel(subm1d), 1);
    s = scatter(subT, subm1d, 64, tcols(idT, :), 'filled');
    
end

npairs = 0;
ndecre = 0;
for idf = 1:numel(allids)
    
    subm1d = am1pdfd(fishids == allids(idf));
    subT = temperatures(fishids == allids(idf));

    if numel(subT) < 2
        continue
    end

    [sortsubT, orii] = sort(subT);
    sortsubm1d = subm1d(orii);

    % Check in how many pairs the mean decreases when temperature increases
    idx = find(sortsubm1d);
    pairs = nchoosek(idx, 2);

    rme = NaN(size(pairs, 1), 1);
    dTe = NaN(size(pairs, 1), 1);
    for k = 1:size(pairs, 1) 
        rme(k) = sortsubm1d(pairs(k, 2))/sortsubm1d(pairs(k, 1));
        dTe(k) = sortsubT(pairs(k, 2)) - sortsubT(pairs(k, 1));
    end
    
    npairs = npairs + size(pairs, 1);
    ndecre = ndecre + sum(rme < 1);

    plot(sortsubT, sortsubm1d, 'LineWidth', 1.5, 'Color', [0.6, 0.6, 0.6]);
end

e = errorbar(allT, Tm1pdfd, Ts1pdfd, 'k', 'LineWidth',2);
e.CapSize = 0;
% s = scatter(allT, Tm1pdfd, 64, tcols, 'filled');
xlabel('temperature (°C)');
ylabel('<m_{L,R}(T)>');
ax = gca;
ax.YGrid = 'off';
ax.XTick = allT;
xlim([17.5, 33.5]);

% --- Save
if ~isempty(fnout)
    save(fnout, 'am1pdfd', 'am1pdfi');
end