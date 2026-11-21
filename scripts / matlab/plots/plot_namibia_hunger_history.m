clear; clc; close all;
paths = project_paths();

baseDir = paths.namibia_results;
indir = fullfile(baseDir,'full');
require_input_files(indir, { ...
    'mean_hunger_timeseries.csv', ...
    'hunger_distribution_snapshots.csv'});
outdir = paths.namibia_figures;
if ~isfolder(outdir), mkdir(outdir); end

fontName = 'Times New Roman';
fontSize = 17;
vh = 0.808976066888251;

%% Figure 1: two time-series panels
meanFile = fullfile(indir,'mean_hunger_timeseries.csv');
if ~isfile(meanFile), error('Missing %s. Supply this CSV separately; the included Python runner does not write it.',meanFile); end
M = readtable(meanFile);
fig = figure('Color','w','Position',[100 100 850 700]);
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax = nexttile;
plot(ax,M.time,M.mean_hunger,'LineWidth',1.8,'Color',[0.12 0.37 0.66]);
yline(ax,vh,'r--','LineWidth',1.25);
ylabel(ax,'Mean hunger level');
title(ax,'Population-weighted mean hunger');
set(ax,'FontName',fontName,'FontSize',fontSize,'XColor','k','YColor','k');
grid(ax,'on'); box(ax,'on'); xlim(ax,[0 max(M.time)]);
ax = nexttile;
plot(ax,M.time,M.fraction_above_threshold,'LineWidth',1.8,'Color',[0.63 0.30 0.12]);
xlabel(ax,'Time (days)'); ylabel(ax,'Fraction above threshold');
title(ax,'Population fraction exceeding the hunger threshold');
set(ax,'FontName',fontName,'FontSize',fontSize,'XColor','k','YColor','k');
grid(ax,'on'); box(ax,'on'); xlim(ax,[0 max(M.time)]); ylim(ax,[0 1]);
exportgraphics(fig,fullfile(outdir,'namibia_mean_hunger_over_time.png'),'Resolution',300);
exportgraphics(fig,fullfile(outdir,'namibia_mean_hunger_over_time.eps'),'ContentType','vector');
savefig(fig,fullfile(outdir,'namibia_mean_hunger_over_time.fig'));

%% Figure 2: five hunger-distribution snapshots
snapFile = fullfile(indir,'hunger_distribution_snapshots.csv');
if ~isfile(snapFile), error('Missing %s. Supply this CSV separately; the included Python runner does not write it.',snapFile); end
S = readtable(snapFile);
days = [0 45 90 135 180];
fig2 = figure('Color','w','Position',[100 100 1600 350]);
tiledlayout(fig2,1,5,'TileSpacing','compact','Padding','compact');
for k = 1:numel(days)
    ax = nexttile;
    col = sprintf('day_%d',days(k));
    if ~ismember(col,S.Properties.VariableNames)
        error('Missing column %s in %s',col,snapFile);
    end
    plot(ax,S.hunger_level,S.(col),'LineWidth',1.8,'Color',[0.12 0.37 0.66]);
    hold(ax,'on'); xline(ax,vh,'r--','LineWidth',1.3);
    xlim(ax,[0 5]); ylim(ax,[0 1]);
    title(ax,sprintf('Day %d',days(k)),'Color', 'k');
    xlabel(ax,'Hunger level');
    if k == 1, ylabel(ax,'Fraction of living population'); end
    set(ax,'FontName',fontName,'FontSize',fontSize,'XColor','k','YColor','k');
    grid(ax,'on'); box(ax,'on');
end
exportgraphics(fig2,fullfile(outdir,'namibia_hunger_distribution_snapshots.png'),'Resolution',300);
exportgraphics(fig2,fullfile(outdir,'namibia_hunger_distribution_snapshots.eps'),'ContentType','vector');
savefig(fig2,fullfile(outdir,'namibia_hunger_distribution_snapshots.fig'));
fprintf('Figures saved to %s\n',outdir);
