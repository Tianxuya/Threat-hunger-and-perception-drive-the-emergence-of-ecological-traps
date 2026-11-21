clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================
baseDir = paths.namibia_results;

indir = fullfile(baseDir, 'full');
require_input_files(indir, {'hunger_distribution_timeseries.csv'});

outdir = paths.namibia_figures;

if ~isfolder(outdir), mkdir(outdir); end

%% ================= LOAD DATA =================
T = readtable(fullfile(indir, ...
    'hunger_distribution_timeseries.csv'), ...
    'VariableNamingRule','preserve');

time = T{:,1};
P = T{:,2:end};   % Rows = time, columns = hunger bins

% Recover hunger-grid values from CSV column names
names = T.Properties.VariableNames(2:end);
v = cellfun(@(x) sscanf(x,'v_%f'), names);

%% ================= SETTINGS =================
fontName = 'Times New Roman';
fontSize = 25;

vh = 0.808976066888251;

% Blue colormap, light to dark
blueAnchors = [
    0.95 0.97 1.00
    0.67 0.82 0.94
    0.31 0.62 0.82
    0.06 0.30 0.58
    ];

cmap = interp1(linspace(0,0.05,4),blueAnchors, ...
    linspace(0,0.05,256));

%% ================= PLOT =================
fig = figure('Color','w','Position',[100 100 1100 500]);
ax = axes(fig);

fig.Color = 'w';
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';

imagesc(ax,time,v,P');

set(ax,'YDir','normal');

colormap(ax,parula(256));
clim(ax,[0 0.05]);

hold(ax,'on');
yline(ax,vh,'r--','LineWidth',1.8);

xlim(ax,[0 180]);
ylim(ax,[0 5]);

xlabel(ax,'Time (days)', ...
    'FontName',fontName,'FontSize',fontSize);

ylabel(ax,'Hunger level', ...
    'FontName',fontName,'FontSize',fontSize);

% title(ax,'Temporal evolution of the hunger distribution', ...
%     'FontName',fontName,'FontSize',fontSize, ...
%     'FontWeight','normal');


%% ================= COLORBAR =================
cb = colorbar(ax);

cb.Ticks = [0 0.01 0.02 0.03 0.04 0.05];
cb.FontName = fontName;
cb.FontSize = 20;

cb.Color = 'k';                % Black tick labels
cb.Label.Color = 'k';          % Black colorbar title
cb.Label.String = 'Fraction of living population';
cb.Label.FontName = fontName;
cb.Label.FontSize = 25;

cb.TickLength = 0;

ax.FontName = fontName;
ax.FontSize = fontSize;
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1;

box(ax,'on');

%% ================= EXPORT =================
exportgraphics(fig,fullfile(outdir, ...
    'namibia_hunger_distribution_heatmap.png'), ...
    'Resolution',300);

exportgraphics(fig,fullfile(outdir, ...
    'namibia_hunger_distribution_heatmap.eps'), ...
    'ContentType','vector');

savefig(fig,fullfile(outdir, ...
    'namibia_hunger_distribution_heatmap.fig'));
