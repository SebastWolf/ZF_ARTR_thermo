%% ================================ SI Figure 2D ============================================
% --- Persistence time vs mean(mL mR)

clear
clc

% --- Definitions
pwd = '/Users/wolf/Documents/CNRS/Papers/Elife 2023'
isinpath = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep 'ARTR_inference_J' filesep];

fn_persis = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep '/Behavior data/persistime.mat'];
fn_pdf00 = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep '/Behavior data/pdf00.mat'];
fn_pdfm1 = [pwd filesep 'Data availability' filesep 'Thermotaxis' filesep '/Behavior data/meanmLmR.mat'];

allT = [18, 22, 26, 30, 33];

% --- Parameters
colors = rainbow(numel(allT));

% --- Load data
fpersis = load(fn_persis);
fpdf00 = load(fn_pdf00);
fpdfm1 = load(fn_pdfm1);

% --- Persistence time vs pdf(0,0)
figure('Name', 'persistimep00'); hold on
[~, I] = ismember(fpersis.temperatures, allT);
cols = colors(I, :);
scatter(fpdf00.p00d, fpersis.mdata, [], cols, 'filled');
xlabel('pdf(0,0)');
ylabel('persistence time (s)');

blank = NaN(numel(allT), 2);
dn = arrayfun(@num2str, allT', 'UniformOutput', false);
s = arrayfun(@(x,y,z1,z2,z3) scatter(x, y, [], [z1, z2, z3], 'filled'), ...
    blank(:, 1), blank(:, 2), colors(:, 1), colors(:, 2), colors(:, 3));

lg = legend(s(:), dn);
title(lg, 'temperature (°C)');
title('Data');

% --- Persistence time vs mean(mL mR)
figure('Name', 'persistimeM1'); hold on
[~, I] = ismember(fpersis.temperatures, allT);
cols = colors(I, :);
scatter(fpdfm1.am1pdfd, fpersis.mdata, [], cols, 'filled');
xlabel('mean(m_L,m_R)');
ylabel('persistence time (s)');

blank = NaN(numel(allT), 2);
dn = arrayfun(@num2str, allT', 'UniformOutput', false);
s = arrayfun(@(x,y,z1,z2,z3) scatter(x, y, [], [z1, z2, z3], 'filled'), ...
    blank(:, 1), blank(:, 2), colors(:, 1), colors(:, 2), colors(:, 3));

ylim([0, 18]);

lg = legend(s(:), dn, 'Location', 'northwest');
title(lg, 'temperature (°C)');
title('Data');
grid off
axis square
ax = gca;
ax.XTick = [0, 0.1, 0.2, 0.3];
R = corrcoef([fpdfm1.am1pdfd', fpersis.mdata]);
text(0.25, 4, ['R = ' num2str(R(1, 2), '%0.2f')], 'FontSize', 12);