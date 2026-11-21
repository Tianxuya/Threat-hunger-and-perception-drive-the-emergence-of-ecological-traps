
%% ============================================================
% BIFURCATION COMPARISON: v0 = 0, 2, 10
%
% Figure 1: Individual patch occupancy (3 panels)
% Figure 2: Total patch occupancy (3 panels)
%
% Shaded regions represent +/- 1 standard deviation.
%
% Original simulation data are unchanged.
%% ============================================================

clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================

indir = paths.agent_results;
require_input_files(indir, { ...
    'bifurcation_all_patches_v0.mat', ...
    'bifurcation_all_patches_v2.mat', ...
    'bifurcation_all_patches_v10.mat'});

outdir = paths.agent_figures;

%% ================= LOAD DATA =================

hunger_values = [0 2 10];

data = cell(3,1);

for h = 1:3

    filename = sprintf( ...
        'bifurcation_all_patches_v%g.mat', hunger_values(h));

    data{h} = load(fullfile(indir, filename));

end

%% ================= SETTINGS =================

colors = [
    0.00 0.00 1.00;   % L-L: blue
    0.85 0.30 0.30;   % M-M: red
    0.00 0.00 0.00    % H-H: black
];

names = {'L-L','M-M','H-H'};

fontName = 'Times New Roman';
fontSize = 22;
lineWidth = 2.5;

% Standard deviation shading
shadeAlpha = 0.20;

% Common y limits
ylim_individual = [0 500];
ylim_total = [0 500];

%% ============================================================
% FIGURE 1: INDIVIDUAL PATCH OCCUPANCY
%% ============================================================

fig1 = figure(1);

set(fig1, ...
    'Visible','on', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 12 4]);

tl1 = tiledlayout(fig1,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

line_handles = gobjects(1,3);

for h = 1:3

    S = data{h};

    R = S.cognition_radii(:);
    M = S.mean_counts;
    SD = S.std_counts;

    ax = nexttile(tl1);
    hold(ax,'on');

    for reg = 1:3

        mean_reg = M(:,reg);
        std_reg = SD(:,reg);

        % Shaded standard deviation
        fill_x = [R; flipud(R)];

        fill_y = [
            mean_reg + std_reg;
            flipud(mean_reg - std_reg)
        ];

        fill(ax, ...
            fill_x, fill_y, colors(reg,:), ...
            'FaceAlpha',shadeAlpha, ...
            'EdgeColor','none', ...
            'HandleVisibility','off');

        % Mean occupancy line
        p = plot(ax, ...
            R, mean_reg, ...
            'Color',colors(reg,:), ...
            'LineWidth',lineWidth);

        if h == 1
            line_handles(reg) = p;
        end

    end

    title(ax, ...
        sprintf('v_0 = %g',hunger_values(h)), ...
        'FontName',fontName, ...
        'FontSize',18, ...
        'FontWeight','normal', ...
        'Interpreter','tex');

    xlabel(ax,'Perception radius', ...
        'FontName',fontName, ...
        'FontSize',fontSize);

    set(ax, ...
        'FontName',fontName, ...
        'FontSize',fontSize, ...
        'LineWidth',1.0, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');

    xlim(ax,[min(R) max(R)]);
    ylim(ax,ylim_individual);

    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridAlpha = 0.15;

    if h > 1
        ax.YTickLabel = [];
    end

    hold(ax,'off');

end

ylabel(tl1, ...
    'Final abundance in patch', ...
    'FontName',fontName, ...
    'FontSize',fontSize+1);

lgd1 = legend(line_handles,names, ...
    'Orientation','horizontal', ...
    'FontName',fontName, ...
    'FontSize',18, ...
    'Box','off');

lgd1.Layout.Tile = 'north';

%% ============================================================
% FIGURE 2: TOTAL OCCUPANCY IN ALL THREE PATCHES
%% ============================================================

fig2 = figure(2);

set(fig2, ...
    'Visible','on', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 12 4]);

tl2 = tiledlayout(fig2,1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for h = 1:3

    S = data{h};

    R = S.cognition_radii(:);

    % Mean total occupancy
    total_mean = sum(S.mean_counts,2);

    %% Correct standard deviation of total occupancy
    %
    % all_counts dimensions:
    % repetitions x regions x cognition radii
    %
    % First sum regions for every repetition.

    if isfield(S,'all_counts')

        A = S.all_counts;

        nRep = size(A,1);
        nR = size(A,3);

        % repetitions x cognition radii
        total_rep = reshape(sum(A,2),nRep,nR);

        total_std = std(total_rep,0,1)';

    else

        error(['The file for v0 = %g does not contain ', ...
               'all_counts. Cannot calculate the exact ', ...
               'standard deviation of total occupancy.'], ...
               hunger_values(h));

    end

    ax = nexttile(tl2);
    hold(ax,'on');

    %% Shaded standard deviation

    fill_x = [R; flipud(R)];

    fill_y = [
        total_mean + total_std;
        flipud(total_mean - total_std)
    ];

    fill(ax, ...
        fill_x,fill_y,[0.7 0.7 0.7], ...
        'FaceAlpha',0.30, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');

    %% Mean total occupancy

    plot(ax, ...
        R,total_mean, ...
        'k-', ...
        'LineWidth',lineWidth, ...
        'Marker','o', ...
        'MarkerSize',5);

    title(ax, ...
        sprintf('v_0 = %g',hunger_values(h)), ...
        'FontName',fontName, ...
        'FontSize',20, ...
        'FontWeight','normal', ...
        'Interpreter','tex');

    xlabel(ax,'Perception radius', ...
        'FontName',fontName, ...
        'FontSize',fontSize);

    set(ax, ...
        'FontName',fontName, ...
        'FontSize',fontSize, ...
        'LineWidth',1.0, ...
        'TickDir','out', ...
        'Box','off', ...
        'Layer','top');

    xlim(ax,[min(R) max(R)]);
    ylim(ax,ylim_total);

    ax.YGrid = 'on';
    ax.XGrid = 'off';
    ax.GridAlpha = 0.15;

    if h > 1
        ax.YTickLabel = [];
    end

    hold(ax,'off');

end

ylabel(tl2, ...
    'Total agents in food patches', ...
    'FontName',fontName, ...
    'FontSize',fontSize-1);

%% ================= SAVE FIGURES =================

if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Figure 1
savefig(fig1,fullfile(outdir, ...
    'bifurcation_individual_patches_comparison.fig'));

exportgraphics(fig1,fullfile(outdir, ...
    'bifurcation_individual_patches_comparison.eps'), ...
    'ContentType','vector');

exportgraphics(fig1,fullfile(outdir, ...
    'bifurcation_individual_patches_comparison.png'), ...
    'Resolution',600);

% Figure 2
savefig(fig2,fullfile(outdir, ...
    'bifurcation_total_patches_comparison.fig'));

exportgraphics(fig2,fullfile(outdir, ...
    'bifurcation_total_patches_comparison.eps'), ...
    'ContentType','vector');

exportgraphics(fig2,fullfile(outdir, ...
    'bifurcation_total_patches_comparison.png'), ...
    'Resolution',600);

close(fig1);
close(fig2);

fprintf('Bifurcation comparison figures saved to:\n%s\n',outdir);

%% ============================================================
% FIGURE 3: TOTAL OCCUPANCY — THREE HUNGER LEVELS
%           IN ONE FIGURE
%% ============================================================

fig3 = figure(3);

set(fig3, ...
    'Visible','on', ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[1 1 7.5 5.5]);

hold on;

% Colors for three hunger levels
hunger_colors = [
    0.00 0.00 1.00;    % v0 = 0: blue
    0.85 0.30 0.30;    % v0 = 2: red
    0.00 0.00 0.00     % v0 = 10: black
];

line_handles3 = gobjects(1,3);

for h = 1:3

    S = data{h};

    R = S.cognition_radii(:);

    % Mean total occupancy
    total_mean = sum(S.mean_counts,2);

    % Exact standard deviation of total occupancy
    A = S.all_counts;

    nRep = size(A,1);
    nR = size(A,3);

    total_rep = reshape(sum(A,2),nRep,nR);

    total_std = std(total_rep,0,1)';

    % Standard deviation shading
    fill_x = [R; flipud(R)];

    fill_y = [
        total_mean + total_std;
        flipud(total_mean - total_std)
    ];

    fill(fill_x,fill_y,hunger_colors(h,:), ...
        'FaceAlpha',0.15, ...
        'EdgeColor','none', ...
        'HandleVisibility','off');

    % Mean occupancy
    % line_handles3(h) = plot( ...
    %     R,total_mean, ...
    %     'Color',hunger_colors(h,:), ...
    %     'LineWidth',lineWidth);
    line_handles3(h) = plot( ...
        R, total_mean, ...
        'Color', hunger_colors(h,:), ...
        'LineWidth', lineWidth, ...
        'LineStyle', '-', ...
        'Marker', 'o', ...
        'MarkerSize', 6);

end

xlabel('Perception radius', ...
    'FontName',fontName, ...
    'FontSize',fontSize);

ylabel('Total agents in food patches', ...
    'FontName',fontName, ...
    'FontSize',fontSize);

set(gca, ...
    'FontName',fontName, ...
    'FontSize',fontSize, ...
    'LineWidth',1.0, ...
    'TickDir','out', ...
    'Box','off', ...
    'Layer','top');

ylim(ylim_total);

grid on;
set(gca,'GridAlpha',0.15);

legend(line_handles3, ...
    {'v_0 = 0','v_0 = 2','v_0 = 10'}, ...
    'Interpreter','tex', ...
    'FontName',fontName, ...
    'FontSize',18, ...
    'Location','best', ...
    'Box','off');

hold off;

%% ================= SAVE FIGURE 3 =================

savefig(fig3,fullfile(outdir, ...
    'bifurcation_total_patches_overlay.fig'));

exportgraphics(fig3,fullfile(outdir, ...
    'bifurcation_total_patches_overlay.eps'), ...
    'ContentType','vector');

exportgraphics(fig3,fullfile(outdir, ...
    'bifurcation_total_patches_overlay.png'), ...
    'Resolution',600);

fprintf('Combined bifurcation figure saved.\n');
