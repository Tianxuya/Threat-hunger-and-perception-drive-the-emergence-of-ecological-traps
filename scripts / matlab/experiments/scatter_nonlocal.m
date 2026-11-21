%% ============================================================
% Simulation summary figure using only all_agent_counts files
% Panel (a): mean number of agents across all perception conditions
% Panel (b): proportions for nonlocal perception only
%% ============================================================

clear; clc; close all;

% Resolve outputs from this file, independent of the current directory.
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(fileparts(script_dir)));
results_dir = getenv('AGENT_RESULTS_DIR');
if isempty(results_dir), results_dir = fullfile(repo_dir, 'results', 'agent'); end
figures_dir = getenv('AGENT_FIGURES_DIR');
if isempty(figures_dir), figures_dir = fullfile(repo_dir, 'figures', 'agent'); end
if ~isfolder(results_dir), mkdir(results_dir); end
if ~isfolder(figures_dir), mkdir(figures_dir); end


%% ---------------- Load data ----------------
data_file = fullfile(results_dir, 'v0_Toxin_agent_counts_50reps.mat');
if ~isfile(data_file)
    error('Missing %s. Run all_conditions_comparison with AGENT_INITIAL_HUNGER set to 0 first.', data_file);
end
loaded = load(data_file, 'all_agent_counts', 'num_agents', 'v_start');
if isfield(loaded, 'v_start') && loaded.v_start ~= 0
    error('File %s contains v_start = %g; expected 0.', data_file, loaded.v_start);
end
v0_data = loaded.all_agent_counts;
v0_num_agents = loaded.num_agents;

data_file = fullfile(results_dir, 'v2_Toxin_agent_counts_50reps.mat');
if ~isfile(data_file)
    error('Missing %s. Run all_conditions_comparison with AGENT_INITIAL_HUNGER set to 2 first.', data_file);
end
loaded = load(data_file, 'all_agent_counts', 'num_agents', 'v_start');
if isfield(loaded, 'v_start') && loaded.v_start ~= 2
    error('File %s contains v_start = %g; expected 2.', data_file, loaded.v_start);
end
v2_data = loaded.all_agent_counts;
v2_num_agents = loaded.num_agents;

data_file = fullfile(results_dir, 'v10_Toxin_agent_counts_50reps.mat');
if ~isfile(data_file)
    error('Missing %s. Run all_conditions_comparison with AGENT_INITIAL_HUNGER set to 10 first.', data_file);
end
loaded = load(data_file, 'all_agent_counts', 'num_agents', 'v_start');
if isfield(loaded, 'v_start') && loaded.v_start ~= 10
    error('File %s contains v_start = %g; expected 10.', data_file, loaded.v_start);
end
v10_data = loaded.all_agent_counts;

num_agents = loaded.num_agents;
if ~isscalar(num_agents) || num_agents < 1 || ...
        v0_num_agents ~= num_agents || v2_num_agents ~= num_agents
    error('All three simulation files must use the same positive num_agents.');
end
if ~isequal(size(v0_data), size(v2_data), size(v10_data)) || ...
        size(v0_data, 1) ~= 3 || size(v0_data, 2) ~= 3
    error('Simulation files must have matching 3-by-3-by-nRep count arrays.');
end

%% ---------------- Settings ----------------
condition_labels = {'Random', 'Local', 'Nonlocal'};
hunger_labels = {'v=0', 'v=2', 'v=10'};

% PNAS-style colors
c_R1  = [0.27 0.46 0.70];   % muted blue
c_R2  = [0.89 0.47 0.20];   % muted orange
c_R3  = [0.17 0.63 0.17];   % muted green
c_out = [0.82 0.82 0.82];   % light gray
bar_colors = [c_R1; c_R2; c_R3; c_out];

region_labels = {'Region 1 (L--L)', ...
                 'Region 2 (M--M)', ...
                 'Region 3 (H--H)', ...
                 'No selection'};

%% ============================================================
% Panel (a): mean number of agents
%% ============================================================

count_v0  = compute_patch_counts(v0_data,  num_agents);
count_v2  = compute_patch_counts(v2_data,  num_agents);
count_v10 = compute_patch_counts(v10_data, num_agents);

plot_data_a = [
    count_v0;
    count_v2;
    count_v10
];

x = [1 2 3, 4.8 5.8 6.8, 8.6 9.6 10.6];
group_centers = [2, 5.8, 9.6];

fig = figure(100); clf;
set(fig, 'Units', 'inches', 'Position', [1, 1, 7.2, 6.2]);
set(fig, 'Color', 'w');

%% ---------------- Panel (a) ----------------
ax1 = subplot(2,1,1);
hold(ax1, 'on');

b = bar(ax1, x, plot_data_a, 'stacked', 'BarWidth', 0.72);

for k = 1:4
    b(k).FaceColor = bar_colors(k,:);
    b(k).EdgeColor = 'none';
end

set(ax1, ...
    'XTick', x, ...
    'XTickLabel', repmat(condition_labels, 1, 3), ...
    'FontName', 'Arial', ...
    'FontSize', 25, ...
    'LineWidth', 0.8, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'Layer', 'top');

xtickangle(ax1, 0);
ylabel(ax1, 'Mean number of agents', 'FontName', 'Arial', 'FontSize', 10);

ylim(ax1, [-80 num_agents]);
yticks(ax1, 0:100:500);

ax1.YGrid = 'on';
ax1.XGrid = 'off';
ax1.GridAlpha = 0.18;

label_y = -55;
for h = 1:3
    text(ax1, group_centers(h), label_y, hunger_labels{h}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end

text(ax1, -0.08, 1.05, '(a)', ...
    'Units', 'normalized', ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

legend(ax1, region_labels, ...
    'Orientation', 'horizontal', ...
    'Location', 'northoutside', ...
    'FontName', 'Arial', ...
    'FontSize', 8.5, ...
    'Box', 'off');

hold(ax1, 'off');

%% ============================================================
% Panel (b): nonlocal perception proportions only
%% ============================================================

cond_nonlocal = 3;

prop_nonlocal = [
    compute_condition_proportions(v0_data,  cond_nonlocal, num_agents);
    compute_condition_proportions(v2_data,  cond_nonlocal, num_agents);
    compute_condition_proportions(v10_data, cond_nonlocal, num_agents)
];

plot_data_b = prop_nonlocal * 100;

ax2 = subplot(2,1,2);
hold(ax2, 'on');

b2 = bar(ax2, 1:3, plot_data_b, 'stacked', 'BarWidth', 0.55);

for k = 1:4
    b2(k).FaceColor = bar_colors(k,:);
    b2(k).EdgeColor = 'none';
end

set(ax2, ...
    'XTick', 1:3, ...
    'XTickLabel', hunger_labels, ...
    'FontName', 'Arial', ...
    'FontSize', 9, ...
    'LineWidth', 0.8, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'Layer', 'top');

ylabel(ax2, 'Patch-selection outcome (%)', ...
    'FontName', 'Arial', ...
    'FontSize', 10);

xlabel(ax2, 'Initial hunger level', ...
    'FontName', 'Arial', ...
    'FontSize', 10);

ylim(ax2, [0 100]);
yticks(ax2, 0:20:100);

ax2.YGrid = 'on';
ax2.XGrid = 'off';
ax2.GridAlpha = 0.18;

title(ax2, 'Nonlocal perception', ...
    'FontName', 'Arial', ...
    'FontSize', 10, ...
    'FontWeight', 'normal');

text(ax2, -0.08, 1.05, '(b)', ...
    'Units', 'normalized', ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

hold(ax2, 'off');

%% ---------------- Export ----------------
exportgraphics(fig, fullfile(figures_dir, 'simulation_summary_counts_and_nonlocal_proportion.pdf'), ...
    'ContentType', 'vector');

exportgraphics(fig, fullfile(figures_dir, 'simulation_summary_counts_and_nonlocal_proportion.png'), ...
    'Resolution', 600);

%% ============================================================
% Helper function 1: counts
%% ============================================================
function count_matrix = compute_patch_counts(agent_counts, num_agents)
    % agent_counts: 3 × 3 × nRep
    % rows: Random, Local, Nonlocal
    % columns: Region 1, Region 2, Region 3

    mean_counts = mean(agent_counts, 3);                 % 3 × 3
    outside_counts = num_agents - sum(mean_counts, 2);   % 3 × 1

    count_matrix = [mean_counts, outside_counts];        % 3 × 4
end

%% ============================================================
% Helper function 2: proportions
%% ============================================================
function prop_vec = compute_condition_proportions(agent_counts, cond_idx, num_agents)
    % output: 1 × 4
    % columns: Region 1, Region 2, Region 3, No selection

    mean_counts = mean(agent_counts(cond_idx, :, :), 3);
    mean_counts = reshape(mean_counts, 1, 3);

    region_props = mean_counts ./ num_agents;
    outside_prop = max(0, 1 - sum(region_props));

    prop_vec = [region_props, outside_prop];
end
