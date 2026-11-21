%% ============================================================
% PNAS-style patch occupancy comparison
%
% Three panels: initial hunger v0 = 0, 2, 10
% Three bars per panel: Random, Local, Nonlocal
% Stacked colors: Region 1, Region 2, Region 3
%
% Error bars: 95% confidence intervals for mean total patch
% occupancy, calculated from the total within EACH repetition.
%
% Save figures without displaying them.
%% ============================================================

clear all; clc; close all;
paths = project_paths();

%% ================= LOAD DATA =================

indir = paths.agent_results;
require_input_files(indir, { ...
    'v0_Toxin_agent_counts_50reps.mat', ...
    'v2_Toxin_agent_counts_50reps.mat', ...
    'v10_Toxin_agent_counts_50reps.mat'});

S0  = load(fullfile(indir, 'v0_Toxin_agent_counts_50reps.mat'));
S2  = load(fullfile(indir, 'v2_Toxin_agent_counts_50reps.mat'));
S10 = load(fullfile(indir, 'v10_Toxin_agent_counts_50reps.mat'));

all_data = {
    S0.all_agent_counts
    S2.all_agent_counts
    S10.all_agent_counts
};

%% ================= SETTINGS =================

num_agents = 500;

condition_labels = {'Random','Local','Nonlocal'};
hunger_values = [0, 2, 10];

% Region colors
c_R1 = [0.27 0.46 0.70];   % Blue
c_R2 = [0.89 0.47 0.20];   % Orange
c_R3 = [0.17 0.63 0.17];   % Green

bar_colors = [c_R1; c_R2; c_R3];

fontName = 'Times New Roman';
fontSize = 22;

% 95% CI; set false to use standard deviation instead
use_CI = true;

%% ================= CALCULATE STATISTICS =================

mean_counts = cell(3,1);
total_means = cell(3,1);
total_errors = cell(3,1);

for h = 1:3

    A = all_data{h};

    % Expected dimensions:
    % conditions x regions x repetitions
    % 3 x 3 x nRep

    assert(ndims(A) == 3, ...
        'Data must be a 3D array.');

    assert(size(A,1) == 3 && size(A,2) == 3, ...
        'Expected 3 conditions and 3 regions.');

    nRep = size(A,3);

    % Mean count in each region
    % Result: 3 conditions x 3 regions
    M = mean(A,3);

    mean_counts{h} = M;

    % Total occupancy inside the three regions:
    % sum regions within each repetition FIRST
    total_rep = reshape(sum(A,2), 3, nRep);

    total_means{h} = mean(total_rep,2);

    % Standard deviation across repetitions
    total_sd = std(total_rep,0,2);

    if use_CI

        % Student-t 95% confidence interval for the mean.
        % Requires Statistics and Machine Learning Toolbox.
        tcrit = tinv(0.975,nRep-1);

        total_errors{h} = ...
            tcrit * total_sd / sqrt(nRep);

    else

        total_errors{h} = total_sd;

    end

end

%% ================= SHARED Y-AXIS LIMIT =================

all_upper = [];

for h = 1:3

    all_upper = [
        all_upper;
        total_means{h} + total_errors{h}
    ];

end

ymax = max(all_upper);

% Round upward to convenient axis limit
ymax = ceil(ymax/50)*50;

%% ================= CREATE FIGURE =================

% Invisible from creation: no figure window is displayed.
fig = figure( ...
    'Visible','on', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 9.5 3.8]);

tl = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

axes_handles = gobjects(3,1);
legend_handles = gobjects(3,1);

%% ================= PLOT THREE PANELS =================

for h = 1:3

    ax = nexttile(tl);
    axes_handles(h) = ax;

    hold(ax,'on');

    M = mean_counts{h};
    x = 1:3;

    % Stacked bars: only the three food-toxin regions
    b = bar(ax,x,M,'stacked', ...
        'BarWidth',0.65);

    for k = 1:3

        b(k).FaceColor = bar_colors(k,:);
        b(k).EdgeColor = 'none';

    end

    % Explicit bar handles keep error bars out of the legend.
    if h == 1
        legend_handles = b(:);
    end

    % Error bars for mean total patch occupancy
    errorbar(ax, ...
        x, ...
        total_means{h}, ...
        total_errors{h}, ...
        'k', ...
        'LineStyle','none', ...
        'LineWidth',1.1, ...
        'CapSize',5);

    % Panel title
    title(ax, ...
        sprintf('v_0 = %g',hunger_values(h)), ...
        'FontName',fontName, ...
        'FontSize',16, ...
        'FontWeight','normal', ...
        'Interpreter','tex');

    % Axes
    set(ax, ...
        'XTick', x, ...
        'XTickLabel', condition_labels, ...
        'XTickLabelRotation', 0, ...  % Horizontal labels
        'FontName', fontName, ...
        'FontSize', fontSize-4, ...
        'LineWidth', 0.8, ...
        'TickDir', 'out', ...
        'Box', 'off', ...
        'Layer', 'top');

    xlim(ax,[0.5 3.5]);
    ylim(ax,[0 ymax]);

    ylabel(ax,'');

    % Grid
    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridAlpha = 0.15;
    ax.GridLineStyle = '-';

    % Shared y axis: labels only on first panel
    if h > 1
        ax.YTickLabel = [];
    end

    hold(ax,'off');

end

%% ================= COMMON Y LABEL =================

ylabel(tl, ...
    'Final abundance', ...
    'FontName',fontName, ...
    'FontSize',fontSize-2);

%% ================= LEGEND =================

lgd = legend(axes_handles(1), legend_handles, ...
    {'Region 1 (L-L)', ...
     'Region 2 (M-M)', ...
     'Region 3 (H-H)'}, ...
    'Orientation','horizontal', ...
    'FontName',fontName, ...
    'FontSize',15, ...
    'Box','off');

lgd.Layout.Tile = 'north';

%% ================= SAVE FIGURE =================

outdir = paths.agent_figures;

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

filename = 'simulation_number_summary';

% Save editable MATLAB figure
savefig(fig, fullfile(outdir, [filename '.fig']));

% Save vector EPS for publication
exportgraphics(fig, ...
    fullfile(outdir, [filename '.eps']), ...
    'ContentType','vector');

% Save high-resolution PNG
exportgraphics(fig, ...
    fullfile(outdir, [filename '.png']), ...
    'Resolution',600);

% close(fig);

fprintf('Figure saved to: %s\n', outdir);
