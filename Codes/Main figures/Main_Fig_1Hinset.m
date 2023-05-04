% Compare kflip from behaviour to nu from ARTR

clear
clc

% --- Definitions
fnkflip = '/home/guillaume/Science/Projects/Neurofish/Programs/ThermoAnalysis/Analysis/Ising/elife/data/kflip.mat';
fnnu = '/home/guillaume/Science/Projects/Neurofish/Programs/ThermoAnalysis/Analysis/Ising/elife/data/artr_nu_psd.mat';

T = [18, 22, 26, 30, 33];
colors = rainbow(numel(T));

% --- Load data
Dkflip = load(fnkflip);
Dnu = load(fnnu);

kflips = Dkflip.kflip_psd;
nus = Dnu.nu;

mkflip = cellfun(@mean, kflips);
ekflip = cellfun(@(x) std(x)./sqrt(numel(x)), kflips);
kneg = mkflip - ekflip;
kpos = mkflip + ekflip;

mnu = cellfun(@(x) mean(x, 'omitnan'), nus);
enu = cellfun(@(x) std(x, 'omitnan')./sqrt(numel(x)), nus);
% mnu = NaN(numel(T), 1);
% enu = NaN(numel(T), 1);
% for idT = 1:numel(T)
%     mnu(idT) = mean(nus(Dnu.temperatures == T(idT)));
%     enu(idT) = std(nus(Dnu.temperatures == T(idT)))./sqrt(sum(Dnu.temperatures == T(idT)));
% end

ft = fit(mkflip, mnu, 'a*x', 'StartPoint', 1);

figure('Name', 'kflipnucomp'); hold on
for idT = 1:numel(T)
    
    e = errorbar(mkflip(idT), mnu(idT), enu(idT), enu(idT), ekflip(idT), ekflip(idT));
    e.CapSize = 0;
    e.LineWidth = 5;
    e.Color = colors(idT, :);
    e.DisplayName = ['T = ' num2str(T(idT)) '°C'];
    e.LineStyle = 'none';
    
end

xlim([0, 0.75]);
ylim([0, 0.75]);
axis square
grid off
xlabel('behavior <k_{flip}> (s^{-1})');
ylabel('ARTR <\nu> (s^{-1})');

plot([0, 0.75], ft.a.*[0, 0.75], 'k--', 'DisplayName', ['slope = ' num2str(ft.a, '%0.2f')]);

legend('Location', 'northwest');