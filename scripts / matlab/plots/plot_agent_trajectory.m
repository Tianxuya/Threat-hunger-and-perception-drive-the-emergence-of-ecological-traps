paths = project_paths();
if ~exist('all_results','var')
    resultFile = fullfile(paths.agent_results, 'food_toxin_max_determined.mat');
    if ~isfile(resultFile), error('Missing simulation result: %s', resultFile); end
    load(resultFile, 'all_results', 'centers', 'num_agents', 'radius_peak', 'x_domain', 'y_domain');
end

%% TRAJECTORY FIGURES (one per condition)
% Assumes: all_results, x_domain, y_domain, centers, radius_peak, num_agents
% Loads simulation results from results/agent when needed.

if ~exist('x_domain','var'), x_domain = -10:0.1:10; end
if ~exist('y_domain','var'), y_domain = -10:0.1:10; end

if ~exist('radius_peak','var'), radius_peak = 3; end
radius = radius_peak;   % circle radius for regions

% Region colours
c_low  = [0.27 0.46 0.70];   % muted blue
c_med  =  [0.89 0.47 0.20]; % muted orange
c_high = [0.17 0.63 0.17];        % muted green

region_colors = {c_low, c_med, c_high};

% Condition colours for trajectories
cond_colors = struct();
cond_colors.Random      = [0.35 0.35 0.35];%[1, 0, 1];   % magenta
cond_colors.NoCognition = [0.35 0.35 0.35];%[0.4, 0.2, 0.7];   % purple
cond_colors.Cognition   = [0.35 0.35 0.35];%[0.0, 0.6, 0.6];   % teal


% Display names for conditions
cond_display = struct();
cond_display.Random      = 'Random';
cond_display.NoCognition = 'Local';
cond_display.Cognition   = 'Nonlocal';

% Which conditions are present in all_results
conds_present = fieldnames(all_results);   % e.g., {'Random','NoCognition','Cognition'}

% -------------------------------------------------------------------------
% 1. Separate figures for each condition
% -------------------------------------------------------------------------
save_dir = paths.agent_figures;
if ~isfolder(save_dir), mkdir(save_dir); end
max_traj_agents = 20;   % number of agents to show trajectories for (will downsample if > num_agents)

for i = 1:numel(conds_present)
    cond = conds_present{i};
    results = all_results.(cond);
    traj_color = cond_colors.(cond);
    display_name = cond_display.(cond);

    fig = figure('Position', [100, 100, 600, 400], 'Color', 'white', 'Visible', 'off');


    % Call the local function (defined at the end)
    plot_condition_trajectories(results, display_name, x_domain, y_domain, centers, radius, ...
        region_colors, traj_color, num_agents, max_traj_agents);

    set(findall(gcf, '-property', 'FontName'), ...
        'FontName', 'Times New Roman');
    % Save
    save_name = fullfile(save_dir, sprintf('%s', lower(cond)));
    saveas(fig, [save_name '.fig']);

    % Save PNG
    exportgraphics( ...
        fig, ...
        [save_name '.png'], ...
        'Resolution', 300);

    % Save EPS
    exportgraphics( ...
        fig, ...
        [save_name '.eps'], ...
        'ContentType', 'vector');

    fprintf('Saved: %s\n', save_name);

end

% Local helper for each condition.
function plot_condition_trajectories(results, display_name, x_domain, y_domain, centers, radius, region_colors, traj_color, num_agents, max_traj)
% results: struct with .positions (num_agents x 2 x num_steps) and .alive
% display_name: string for title
% traj_color: RGB triplet for trajectory lines
% max_traj: number of agents to plot trajectories for (e.g., 20)
hold on;
% Light background
rectangle('Position', [x_domain(1), y_domain(1), diff(x_domain([1,end])), diff(y_domain([1,end]))], ...
    'FaceColor', [0.97,0.97,0.97], 'EdgeColor', 'none');

% Draw semi‑transparent region circles

% Plot a subset of trajectories
step = max(1, floor(num_agents / max_traj));
for i = 1:step:num_agents
    if results.alive(i,end)
        traj = squeeze(results.positions(i,:,:))';
        if size(traj,1) > 1
            plot(traj(:,1), traj(:,2), '-', 'LineWidth', 1.5, 'Color', traj_color);
        end
    end
end

% Collect all start and end points
start_x = []; start_y = []; end_x = []; end_y = [];
for i = 1:num_agents
    if results.alive(i,end)
        traj = squeeze(results.positions(i,:,:))';
        if size(traj,1) > 1
            start_x(end+1) = traj(1,1); start_y(end+1) = traj(1,2);
            end_x(end+1)   = traj(end,1); end_y(end+1)   = traj(end,2);
        end
    end
end
% Plot all start (green circles) and end (red squares)
scatter(start_x, start_y, 20, 'o', 'MarkerEdgeColor', [0.4 0.2 0.7], ...
    'MarkerFaceColor', [0.4 0.2 0.7], 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
scatter(end_x, end_y, 60, 's', 'MarkerEdgeColor', [0.5,0,0], ...
        'MarkerFaceColor', [0.95,0,0], 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);

theta = linspace(0,2*pi,100);
for k = 1:3
    cx = centers{k}(1); cy = centers{k}(2);
    circ_x = cx + radius*cos(theta);
    circ_y = cy + radius*sin(theta);
    fill(circ_x, circ_y, region_colors{k}, 'EdgeColor', region_colors{k}, ...
        'FaceAlpha', 0.1, 'LineWidth', 3);
end


% Mark region centres with triangles (using region colours)
for k = 1:3
    plot(centers{k}(1), centers{k}(2), '^', 'MarkerSize', 24, 'LineWidth', 3, ...
        'MarkerFaceColor', region_colors{k}, 'MarkerEdgeColor', 'k');
end
% Text labels near centres
labels = {'L-L', 'M-M', 'H-H'};


for k = 1:3
    text(centers{k}(1), centers{k}(2)-1.5, labels{k}, ...
        'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'BackgroundColor', [1,1,1], 'EdgeColor', region_colors{k}, 'LineWidth', 1.5);
end

xlabel('$x$', 'FontSize', 40, 'FontWeight', 'bold',Interpreter='latex');
ylabel('$y$', 'FontSize', 45, 'FontWeight', 'bold',Interpreter='latex');
% title(display_name, 'FontSize', 32, 'FontWeight', 'bold');
xlim([x_domain(1), x_domain(end)]); ylim([y_domain(1), y_domain(end)]);
axis equal; grid on; box on;
set(gca, 'FontSize', 24, 'LineWidth', 1);

% Legend
%% Point colors
start_color = [0.4 0.2 0.7];   % purple
end_color   = [0.85 0.10 0.10]; % red
h1 = plot(NaN,NaN,'-','LineWidth',3,'Color',traj_color);
h2 = scatter(NaN,NaN,60,'o','MarkerEdgeColor',start_color,'MarkerFaceColor',start_color);
h3 = scatter(NaN,NaN,60,'s','MarkerEdgeColor',end_color,'MarkerFaceColor',end_color);
h4 = plot(NaN,NaN,'^','MarkerSize',18,'LineWidth',2, ...
    'MarkerFaceColor',region_colors{1},'MarkerEdgeColor','k');
h5 = plot(NaN,NaN,'^','MarkerSize',18,'LineWidth',2, ...
    'MarkerFaceColor',region_colors{2},'MarkerEdgeColor','k');
h6 = plot(NaN,NaN,'^','MarkerSize',18,'LineWidth',2, ...
    'MarkerFaceColor',region_colors{3},'MarkerEdgeColor','k');

legend([h1 h2 h3 h4 h5 h6], ...
    {'Agent path','Start points','End points','Region 1','Region 2','Region 3'}, ...
    'Location','best','FontSize',18,'Box','on');
% legend([h1,h2,h4,h5,h6], ...
%     {'Agent path','Start points','Region 1','Region 2','Region 3'}, ...
%     'Location','best', 'FontSize', 18, 'Box','on');

hold off;
end

