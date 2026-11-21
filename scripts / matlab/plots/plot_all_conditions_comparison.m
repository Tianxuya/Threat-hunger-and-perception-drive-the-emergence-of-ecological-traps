paths = project_paths();
require_input_files(paths.agent_results, { ...
    'v0_Toxin_agent_counts_50reps.mat', ...
    'v2_Toxin_agent_counts_50reps.mat', ...
    'v2_noToxin_agent_counts_50reps.mat', ...
    'v10_Toxin_agent_counts_50reps.mat'});

%% BOX PLOT: min, Q1, median, Q3, max for each condition & region
% Load pre‑computed data (after running 50 repetitions once)
% load('agent_counts_50reps.mat');   % loads all_agent_counts, conditions, num_agents
  % load('v2_Toxin_agent_counts_50reps.mat');   % loads all_agent_counts, conditions, num_agents
% load('v2_noToxin_agent_counts_50reps.mat');   % loads all_agent_counts, conditions, num_agents
 % load('v0_Toxin_agent_counts_50reps.mat');   % loads all_agent_counts, conditions, num_agents
% load('v10_noToxin_agent_counts_50reps.mat');   % loads all_agent_counts, conditions, num_agents
load(fullfile(paths.agent_results, 'v0_Toxin_agent_counts_50reps.mat'));

% Colors for regions (matching earlier scheme)
% c_low  = [0 0 1];      % blue
% c_med  = [0.85 0.3 0.3]; % orange
% c_high = [0 0 0];      % red

c_low  = [0.27 0.46 0.70];
c_med  = [0.89 0.47 0.20];
c_high = [0.17 0.63 0.17];
c_out  = [0.82 0.82 0.82];

colors = [c_low; c_med; c_high];

% Figure settings
figure(100);
set(gcf, 'Position', [100, 100, 1200, 600]);
set(gcf, 'DefaultAxesFontSize', 25);
hold on;

% Group positions and offsets
conditions_list = conditions;         % {'Random','NoCognition','Cognition'}
group_positions = 1:length(conditions_list);
bar_width = 0.25;
offsets = [-bar_width, 0, bar_width];

box_handles = gobjects(1,3);

for reg = 1:3
    % Extract data: 50 values for each condition (rows = reps, columns = conditions)
    data_matrix = squeeze(all_agent_counts(:, reg, :))';  % 50 × 3
    
    for cond = 1:3
        data = data_matrix(:, cond);
        data = data(~isnan(data));
        if isempty(data); continue; end
        
        % Compute statistics
        min_val = min(data);
        q1 = quantile(data, 0.25);
        med = median(data);
        q3 = quantile(data, 0.75);
        max_val = max(data);
        
        % x position for this box
        x_center = group_positions(cond) + offsets(reg);
        box_width = bar_width * 0.8;
        
        % Draw box from Q1 to Q3
        rectangle('Position', [x_center - box_width/2, q1, box_width, q3 - q1], ...
                  'FaceColor', colors(reg,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
        
        % Median line
        line([x_center - box_width/2, x_center + box_width/2], [med, med], ...
             'Color', 'k', 'LineWidth', 2.5);
        
        % Whiskers (min and max)
        line([x_center, x_center], [min_val, q1], 'Color', 'k', 'LineWidth', 1.5);
        line([x_center, x_center], [q3, max_val], 'Color', 'k', 'LineWidth', 1.5);
        
        % Caps at whisker ends
        cap_half = box_width/4;
        line([x_center - cap_half, x_center + cap_half], [min_val, min_val], ...
             'Color', 'k', 'LineWidth', 1.5);
        line([x_center - cap_half, x_center + cap_half], [max_val, max_val], ...
             'Color', 'k', 'LineWidth', 1.5);
        
        % Optional: mark mean with a small circle
        mean_val = mean(data);
        plot(x_center, mean_val, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'w');
    end
    
    % Dummy patch for legend
    box_handles(reg) = patch(NaN, NaN, colors(reg,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
end

% Axes and labels
set(gca, 'XTick', group_positions, 'XTickLabel', {'Random', 'Local perception', 'Nonlocal perception'});
ylabel('Number of agents');
legend(box_handles, {'Region 1 (low-low)', 'Region 2 (med-med)', 'Region 3 (high-high)'}, ...
       'Location', 'northwest', 'FontSize', 25);
grid on; box on;
hold off;


%% FIGURE 5: Pie charts using mean of 50 repetitions (legend at bottom, shifted left)

% Compute mean counts across repetitions (3×3 matrix)
mean_counts = mean(all_agent_counts, 3);

figure(5);
clf;
set(gcf, 'Position', [100, 100, 1400, 800]);   % tall figure

% Default font sizes
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultTextFontSize', 25);

% Compute proportions from mean counts
region_prop = mean_counts ./ num_agents;      % 3×3
outside_prop = 1 - sum(region_prop, 2);       % 3×1

% Colors (same as before)
% c_low  = [0 0 1];        % blue
% c_med  = [0.85 0.3 0.3]; % orange
% c_high = [0 0 0];        % red
% c_out  = [0.5 0.5 0.5];  % gray
c_low  = [0.27 0.46 0.70];
c_med  = [0.89 0.47 0.20];
c_high = [0.17 0.63 0.17];
c_out  = [0.82 0.82 0.82];

pie_colors = [c_low; c_med; c_high; c_out];

% Display names for conditions (define if not already in workspace)
display_names = {'Random', 'No cognition', 'Cognition'};

for cond = 1:3
    subplot(1, 3, cond);
    
    data = [region_prop(cond,1), region_prop(cond,2), region_prop(cond,3), outside_prop(cond)];
    percent = data * 100;
    
    % Labels: only percentages
    labels = {sprintf('%.1f%%', percent(1)), ...
              sprintf('%.1f%%', percent(2)), ...
              sprintf('%.1f%%', percent(3)), ...
              sprintf('%.1f%%', percent(4))};
    
    % Pie chart
    p = pie(data, labels);
    colormap(pie_colors);
    
    % Title moved down (negative y position)
    th = title(display_names{cond}, 'FontSize', 25);
    th.Position(2) = -1.8;
    
    % Set pie label font size (22 as in your snippet)
    text_objs = findobj(p, 'Type', 'text');
    set(text_objs, 'FontSize', 22);
end

% Legend at bottom, horizontal, shifted left
lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
             'Location', 'southoutside', 'Orientation', 'horizontal', ...
             'FontSize', 25, 'Box', 'off');

% Shift legend left and down (as you had)
lgd.Position(1) = lgd.Position(1) - 0.15;
lgd.Position(2) = lgd.Position(2) - 0.15;

drawnow;

%% PIE CHARTS: Compare no toxin vs with toxin for specific conditions
% Load data (both files must exist)
load(fullfile(paths.agent_results, 'v2_noToxin_agent_counts_50reps.mat'));   % variable: all_agent_counts
noToxin_data = all_agent_counts;           % 3 conditions × 3 regions × 50 reps

load(fullfile(paths.agent_results, 'v2_Toxin_agent_counts_50reps.mat')); % variable: all_agent_counts
withToxin_data = all_agent_counts;

% Define condition indices: 1=Random, 2=NoCognition, 3=Cognition
cond_cog = 3;
cond_nocog = 2;

% Number of agents (should be same in both files)
num_agents = 500;   % adjust if different

% Colors (same scheme)
c_low  = [0 0 1];       % blue
c_med  = [0.85 0.3 0.3];% orange
c_high = [0 0 0];       % red
c_out  = [0.5 0.5 0.5]; % gray
pie_colors = [c_low; c_med; c_high; c_out];

% Display names for the two toxin conditions
toxin_labels = {'No toxin', 'With toxin'};

%% FIGURE 6: Cognition – no toxin vs with toxin (pie charts)
figure(6);
clf;
set(gcf, 'Position', [100, 100, 1200, 600]);
set(gcf, 'DefaultAxesFontSize', 25);

% Compute mean counts for Cognition
mean_noToxin_cog = mean(noToxin_data(cond_cog, :, :), 3);   % 1×3 vector
mean_withToxin_cog = mean(withToxin_data(cond_cog, :, :), 3);

% Proportions
prop_noToxin = mean_noToxin_cog / num_agents;
prop_withToxin = mean_withToxin_cog / num_agents;
outside_noToxin = 1 - sum(prop_noToxin);
outside_withToxin = 1 - sum(prop_withToxin);

% Data for each pie: [R1, R2, R3, outside]
data_noToxin = [prop_noToxin, outside_noToxin];
data_withToxin = [prop_withToxin, outside_withToxin];

% Create subplots
for idx = 1:2
    subplot(1, 2, idx);
    if idx == 1
        data = data_noToxin;
        title_str = toxin_labels{1};
    else
        data = data_withToxin;
        title_str = toxin_labels{2};
    end
    percent = data * 100;
    labels = {sprintf('%.1f%%', percent(1)), ...
              sprintf('%.1f%%', percent(2)), ...
              sprintf('%.1f%%', percent(3)), ...
              sprintf('%.1f%%', percent(4))};
    p = pie(data, labels);
    colormap(pie_colors);
    % Title moved down (negative y position)
    th = title(title_str, 'FontSize', 25);
    th.Position(2) = -1.5;
    
    % Set pie label font size (22 as in your snippet)
    text_objs = findobj(p, 'Type', 'text');
    set(text_objs, 'FontSize', 22);
end

% Legend at bottom, horizontal, shifted left
lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
             'Location', 'southoutside', 'Orientation', 'horizontal', ...
             'FontSize', 25, 'Box', 'off');

% Shift legend left and down (as you had)
lgd.Position(1) = lgd.Position(1)-0.12 ;
lgd.Position(2) = lgd.Position(2) -0.08;
% Overall title
% sgtitle('Cognition condition: effect of toxin', 'FontSize', 32, 'FontWeight', 'bold');


%% FIGURE 7: NoCognition – no toxin vs with toxin (pie charts)
figure(7);
clf;
set(gcf, 'Position', [100, 100, 1200, 600]);
set(gcf, 'DefaultAxesFontSize', 25);

mean_noToxin_nocog = mean(noToxin_data(cond_nocog, :, :), 3);
mean_withToxin_nocog = mean(withToxin_data(cond_nocog, :, :), 3);

prop_noToxin = mean_noToxin_nocog / num_agents;
prop_withToxin = mean_withToxin_nocog / num_agents;
outside_noToxin = 1 - sum(prop_noToxin);
outside_withToxin = 1 - sum(prop_withToxin);

data_noToxin = [prop_noToxin, outside_noToxin];
data_withToxin = [prop_withToxin, outside_withToxin];

for idx = 1:2
    subplot(1, 2, idx);
    if idx == 1
        data = data_noToxin;
        title_str = toxin_labels{1};
    else
        data = data_withToxin;
        title_str = toxin_labels{2};
    end
    percent = data * 100;
    labels = {sprintf('%.1f%%', percent(1)), ...
              sprintf('%.1f%%', percent(2)), ...
              sprintf('%.1f%%', percent(3)), ...
              sprintf('%.1f%%', percent(4))};
    p = pie(data, labels);
    colormap(pie_colors);
    % Title moved down (negative y position)
    th = title(title_str, 'FontSize', 25);
    th.Position(2) = -1.5;
    
    % Set pie label font size (22 as in your snippet)
    text_objs = findobj(p, 'Type', 'text');
    set(text_objs, 'FontSize', 22);
end

% Legend at bottom, horizontal, shifted left
lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
             'Location', 'southoutside', 'Orientation', 'horizontal', ...
             'FontSize', 25, 'Box', 'off');

% Shift legend left and down (as you had)
lgd.Position(1) = lgd.Position(1)-0.12 ;
lgd.Position(2) = lgd.Position(2) -0.08;

%% PIE CHARTS: Compare low hunger (v=2) vs high hunger (v=10) for specific conditions
load(fullfile(paths.agent_results, 'v2_Toxin_agent_counts_50reps.mat')); % low hunger (v=2)
lowHunger_data = all_agent_counts;       % 3×3×50

load(fullfile(paths.agent_results, 'v10_Toxin_agent_counts_50reps.mat')); % high hunger (v=10) – with toxin
highHunger_data = all_agent_counts;

% Define condition indices: 1=Random, 2=NoCognition, 3=Cognition
cond_cog = 3;
cond_nocog = 2;

% Number of agents (should be same in both files)
num_agents = 500;

% Colors (same scheme)
c_low  = [0 0 1];       % blue
c_med  = [0.85 0.3 0.3]; % orange
c_high = [0 0 0];       % red
c_out  = [0.5 0.5 0.5]; % gray
pie_colors = [c_low; c_med; c_high; c_out];

% Labels for the two hunger conditions
hunger_labels = {'Low hunger', 'High hunger'};

%% FIGURE 8: Cognition – low hunger vs high hunger
figure(8);
clf;
set(gcf, 'Position', [100, 100, 1200, 600]);
set(gcf, 'DefaultAxesFontSize', 25);

% Mean counts for Cognition
mean_low = mean(lowHunger_data(cond_cog, :, :), 3);   % 1×3 vector
mean_high = mean(highHunger_data(cond_cog, :, :), 3);

% Proportions
prop_low = mean_low / num_agents;
prop_high = mean_high / num_agents;
outside_low = 1 - sum(prop_low);
outside_high = 1 - sum(prop_high);

data_low = [prop_low, outside_low];
data_high = [prop_high, outside_high];

% Create subplots
for idx = 1:2
    subplot(1, 2, idx);
    if idx == 1
        data = data_low;
        title_str = hunger_labels{1};
    else
        data = data_high;
        title_str = hunger_labels{2};
    end
    percent = data * 100;
    labels = {sprintf('%.1f%%', percent(1)), ...
              sprintf('%.1f%%', percent(2)), ...
              sprintf('%.1f%%', percent(3)), ...
              sprintf('%.1f%%', percent(4))};
    p = pie(data, labels);
    colormap(pie_colors);
    th = title(title_str, 'FontSize', 25);
    th.Position(2) = -1.5;
    text_objs = findobj(p, 'Type', 'text');
    set(text_objs, 'FontSize', 22);
end

% Legend at bottom
lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
             'Location', 'southoutside', 'Orientation', 'horizontal', ...
             'FontSize', 25, 'Box', 'off');
lgd.Position(1) = lgd.Position(1) - 0.12;
lgd.Position(2) = lgd.Position(2) - 0.08;

% Overall subtitle (optional)
% sgtitle('Cognition condition: effect of initial hunger', 'FontSize', 28, 'FontWeight', 'bold');
% 
%% FIGURE 9: NoCognition – low hunger vs high hunger
figure(9);
clf;
set(gcf, 'Position', [100, 100, 1200, 600]);
set(gcf, 'DefaultAxesFontSize', 25);

mean_low = mean(lowHunger_data(cond_nocog, :, :), 3);
mean_high = mean(highHunger_data(cond_nocog, :, :), 3);

prop_low = mean_low / num_agents;
prop_high = mean_high / num_agents;
outside_low = 1 - sum(prop_low);
outside_high = 1 - sum(prop_high);

data_low = [prop_low, outside_low];
data_high = [prop_high, outside_high];

for idx = 1:2
    subplot(1, 2, idx);
    if idx == 1
        data = data_low;
        title_str = hunger_labels{1};
    else
        data = data_high;
        title_str = hunger_labels{2};
    end
    percent = data * 100;
    labels = {sprintf('%.1f%%', percent(1)), ...
              sprintf('%.1f%%', percent(2)), ...
              sprintf('%.1f%%', percent(3)), ...
              sprintf('%.1f%%', percent(4))};
    p = pie(data, labels);
    colormap(pie_colors);
    th = title(title_str, 'FontSize', 25);
    th.Position(2) = -1.5;
    text_objs = findobj(p, 'Type', 'text');
    set(text_objs, 'FontSize', 22);
end

lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
             'Location', 'southoutside', 'Orientation', 'horizontal', ...
             'FontSize', 25, 'Box', 'off');
lgd.Position(1) = lgd.Position(1) - 0.12;
lgd.Position(2) = lgd.Position(2) - 0.08;

% sgtitle('No cognition condition: effect of initial hunger', 'FontSize', 28, 'FontWeight', 'bold');
%% COMPARE ALL THREE HUNGER LEVELS IN ONE FIGURE (Cognition & NoCognition)
% Load data files
load(fullfile(paths.agent_results, 'v0_Toxin_agent_counts_50reps.mat'));   % fully fed (v=0)
v0_data = all_agent_counts;                % 3×3×50

load(fullfile(paths.agent_results, 'v2_Toxin_agent_counts_50reps.mat'));   % low hunger (v=2)
v2_data = all_agent_counts;

load(fullfile(paths.agent_results, 'v10_Toxin_agent_counts_50reps.mat'));  % high hunger (v=10)
v10_data = all_agent_counts;

num_agents = 500;
cond_cog = 3;      % Cognition
cond_nocog = 2;    % NoCognition

% Colors (same as before)
c_low  = [0 0 1];        % blue
c_med  = [0.85 0.3 0.3]; % orange
c_high = [0 0 0];        % red
c_out  = [0.5 0.5 0.5];  % gray
pie_colors = [c_low; c_med; c_high; c_out];

% --- Figure 8: Cognition: v0 vs v2 vs v10 ---
create_triple_pie_figure(v0_data(cond_cog,:,:), v2_data(cond_cog,:,:), v10_data(cond_cog,:,:), ...
                         'Cognition: effect of initial hunger', num_agents, pie_colors);

% --- Figure 9: NoCognition: v0 vs v2 vs v10 ---
create_triple_pie_figure(v0_data(cond_nocog,:,:), v2_data(cond_nocog,:,:), v10_data(cond_nocog,:,:), ...
                         'No Cognition: effect of initial hunger', num_agents, pie_colors);

% Helper: create a 3‑pie figure for a given condition
function create_triple_pie_figure(data_v0, data_v2, data_v10, condition_name, num_agents, colors)
    figure();
    set(gcf, 'Position', [100, 100, 1600, 600]);   % wider for three pies
    set(gcf, 'DefaultAxesFontSize', 25);
    
    % Compute mean counts and proportions for each hunger level
    mean_v0 = mean(data_v0, 3);     % 1×3
    mean_v2 = mean(data_v2, 3);
    mean_v10 = mean(data_v10, 3);
    
    prop_v0 = mean_v0 / num_agents;
    prop_v2 = mean_v2 / num_agents;
    prop_v10 = mean_v10 / num_agents;
    
    outside_v0 = 1 - sum(prop_v0);
    outside_v2 = 1 - sum(prop_v2);
    outside_v10 = 1 - sum(prop_v10);
    
    data_v0_pie = [prop_v0, outside_v0];
    data_v2_pie = [prop_v2, outside_v2];
    data_v10_pie = [prop_v10, outside_v10];
    
    hunger_labels = {'Fully fed', 'Low hunger', 'High hunger'};
    all_data = {data_v0_pie, data_v2_pie, data_v10_pie};
    
    for idx = 1:3
        subplot(1, 3, idx);
        data = all_data{idx};
        percent = data * 100;
        labels = {sprintf('%.1f%%', percent(1)), ...
                  sprintf('%.1f%%', percent(2)), ...
                  sprintf('%.1f%%', percent(3)), ...
                  sprintf('%.1f%%', percent(4))};
        p = pie(data, labels);
        colormap(colors);
        th = title(hunger_labels{idx}, 'FontSize', 25);
        th.Position(2) = -1.8;
        
        text_objs = findobj(p, 'Type', 'text');
        set(text_objs, 'FontSize', 22);
    end
    
    % Shared legend at bottom
    lgd = legend({'Region 1', 'Region 2', 'Region 3', 'No selection'}, ...
                 'Location', 'southoutside', 'Orientation', 'horizontal', ...
                 'FontSize', 25, 'Box', 'off');
    lgd.Position(1) = lgd.Position(1) - 0.12;
    lgd.Position(2) = lgd.Position(2) - 0.06;
    % sgtitle(condition_name, 'FontSize', 32, 'FontWeight', 'bold');
end

