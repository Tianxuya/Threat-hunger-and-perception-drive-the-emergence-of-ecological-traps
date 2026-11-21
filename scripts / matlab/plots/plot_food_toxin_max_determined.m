paths = project_paths();
if ~exist('all_results','var')
    resultFile = fullfile(paths.agent_results, 'food_toxin_max_determined.mat');
    if ~isfile(resultFile), error('Missing simulation result: %s', resultFile); end
    load(resultFile, 'all_results', 'conditions', 'centers', 'num_agents');
end

%% ============================================================
% Final agent positions: original two-figure plotting workflow
%
% Run food_toxin_max_determined.m BEFORE this plotting script.
%
% Required workspace variables:
%   all_results, conditions, centers, num_agents
%
% Do NOT clear the workspace here.
%
% Font: Times New Roman
% Agent dots: dark gray
%% ============================================================

%% ================= CHECK WORKSPACE =================

if ~exist('all_results','var') || ~exist('conditions','var') || ...
        ~exist('centers','var') || ~exist('num_agents','var')

    error(['Run food_toxin_max_determined.m first, or load its saved ', ...
        'workspace variables: all_results, conditions, centers, num_agents.']);

end

%% ================= FONT AND DOT COLOR =================

fontName = 'Times New Roman';

dot_color = [0.25 0.25 0.25];   % Dark gray

%% ================= ORIGINAL SETTINGS =================

x_min = -10; x_max = 10; dx = 0.1;
y_min = -10; y_max = 10; dy = 0.1;

x_domain = x_min:dx:x_max;
y_domain = y_min:dy:y_max;

% Original region colors
c_low  = [0.27 0.46 0.70];
c_med  = [0.89 0.47 0.20];
c_high = [0.17 0.63 0.17];
c_out  = [0.82 0.82 0.82];

pie_colors = [c_low; c_med; c_high; c_out];

x_lim = [x_domain(1), x_domain(end)];
y_lim = [y_domain(1), y_domain(end)];

display_names = {'Random', 'No cognition', 'Cognition'};

%% ================= FIGURE 1: ALL CONDITIONS =================
% 
fig1 = figure(1);
clf(fig1);

set(fig1, ...
    'Position', [100, 100, 1400, 400], ...
    'DefaultAxesFontSize', 25, ...
    'DefaultAxesFontName', fontName, ...
    'DefaultTextFontName', fontName);

for cond_idx = 1:length(conditions)

    condition = conditions{cond_idx};
    results = all_results.(condition);

    subplot(1, 3, cond_idx);
    hold on;

    %% Collect final positions of living agents

    end_points_x = [];
    end_points_y = [];

    for i = 1:num_agents

        if results.alive(i, end)

            final_x = results.positions(i, 1, end);
            final_y = results.positions(i, 2, end);

            end_points_x(end+1) = final_x;
            end_points_y(end+1) = final_y;

        end

    end

    if ~isempty(end_points_x)

        %% Agent dots
        % Original size 45 and transparency 0.6 are unchanged.

        scatter(end_points_x, end_points_y, ...
            45, dot_color, ...
            'filled', ...
            'MarkerFaceAlpha', 0.6);

        %% Original patch-center stars

        plot(centers{1}(1), centers{1}(2), '*', ...
            'Color', c_low, ...
            'MarkerSize', 34, ...
            'LineWidth', 3);

        plot(centers{2}(1), centers{2}(2), '*', ...
            'Color', c_med, ...
            'MarkerSize', 34, ...
            'LineWidth', 3);

        plot(centers{3}(1), centers{3}(2), '*', ...
            'Color', c_high, ...
            'MarkerSize', 34, ...
            'LineWidth', 3);

        %% Original formatting, with updated font

        xlabel('x', ...
            'FontName', fontName, ...
            'FontSize', 25);

        ylabel('y', ...
            'FontName', fontName, ...
            'FontSize', 25);

        title(display_names{cond_idx}, ...
            'FontName', fontName, ...
            'FontSize', 28, ...
            'FontWeight', 'bold');

        xlim(x_lim);
        ylim(y_lim);

        axis equal;
        grid on;

        set(gca, ...
            'FontName', fontName, ...
            'FontSize', 22, ...
            'LineWidth', 2);

    else

        text(0, 0, 'No agents alive', ...
            'HorizontalAlignment', 'center', ...
            'FontName', fontName, ...
            'FontSize', 25);

        xlim(x_lim);
        ylim(y_lim);

        title(display_names{cond_idx}, ...
            'FontName', fontName, ...
            'FontSize', 28, ...
            'FontWeight', 'bold');

    end

    hold off;

end

%% ================= FIGURE 2: SELECTED CONDITION =================

fig2 = figure(2);
clf(fig2);

set(fig2, ...
    'Position', [100, 100, 650, 600], ...
    'DefaultAxesFontSize', 25, ...
    'DefaultAxesFontName', fontName, ...
    'DefaultTextFontName', fontName);

% Retain your original data selection.
% 'NoCognition' = Local
% 'Cognition'   = Nonlocal
% 'Random'      = Random

condition = 'Cognition';

results = all_results.(condition);

% Match the title to the selected condition.
switch condition

    case 'Random'
        single_title = 'Random';

    case 'NoCognition'
        single_title = 'Local';

    case 'Cognition'
        single_title = 'Nonlocal';

    otherwise
        single_title = condition;

end

hold on;

%% Collect final positions of living agents

end_points_x = [];
end_points_y = [];

for i = 1:num_agents

    if results.alive(i, end)

        final_x = results.positions(i, 1, end);
        final_y = results.positions(i, 2, end);

        end_points_x(end+1) = final_x;
        end_points_y(end+1) = final_y;

    end

end

%% Plot final positions

if ~isempty(end_points_x)

    scatter(end_points_x, end_points_y, ...
        45, dot_color, ...
        'filled', ...
        'MarkerFaceAlpha', 0.9);

    %% Original patch-center stars

    plot(centers{1}(1), centers{1}(2), '*', ...
        'Color', c_low, ...
        'MarkerSize', 34, ...
        'LineWidth', 3);

    plot(centers{2}(1), centers{2}(2), '*', ...
        'Color', c_med, ...
        'MarkerSize', 34, ...
        'LineWidth', 3);

    plot(centers{3}(1), centers{3}(2), '*', ...
        'Color', c_high, ...
        'MarkerSize', 34, ...
        'LineWidth', 3);

else

    text(mean(x_lim), mean(y_lim), ...
        'No agents alive', ...
        'HorizontalAlignment', 'center', ...
        'FontName', fontName, ...
        'FontSize', 40);

end

%% Original formatting, with updated font

% xlabel('x', ...
%     'FontName', fontName, ...
%     'FontSize', 25);
% 
% ylabel('y', ...
%     'FontName', fontName, ...
%     'FontSize', 25);

% title(single_title, ...
%     'FontName', fontName, ...
%     'FontSize', 28, ...
%     'FontWeight', 'bold');

xlim(x_lim);
ylim(y_lim);

axis equal;
grid on;
box on;

set(gca, ...
    'FontName', fontName, ...
    'FontSize', 40, ...
    'LineWidth', 2);



outdir = paths.agent_figures;
if ~isfolder(outdir), mkdir(outdir); end

saveas(fig2, fullfile(outdir, 'scatter_v2_nonlocal.eps'), 'epsc');

hold off;