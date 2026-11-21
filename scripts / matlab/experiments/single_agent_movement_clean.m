% ONE AGENT: LINK MOVEMENT TO THE HUNGER AND TOXIN DECISION TERMS
% Run this script directly. It implements the original Cognition condition,
% with one agent, v(0)=v_start, and the fixed initial position (0,0).
%
% Optional data/agent/single_agent_movement_input.csv (from the repository
% root) supplies standard-normal draws. Its repeated initial_x/initial_y
% values must both be zero. The
% Cognition trajectory is simulated anew from the centre. If the CSV is
% absent, rng(random_seed,'twister') supplies reproducible normal draws; the initial
% position remains (0,0). MATLAB and other languages need not share an RNG.
%
% Timing matches the source: update hunger using food at the OLD position;
% then compare H_n=d_v*v_{n+1} and T_n=d_tau*tau(x_n,y_n).
% The figure uses ACTUAL nonzero directed movement to colour each interval.
% An eligible food/avoidance rule can still produce diffusion only when the
% perceived field is uniform or the selected target is the current position.
% Both figures explicitly use a light theme and white export backgrounds,
% independently of MATLAB's desktop theme or the operating system theme.

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


%% Parameters and output location
use_saved_random_input = true;
random_seed = 6;
num_steps = 500; dt = 0.05;        % 500 states; 499 decisions; final t=24.95
D = 1; p_s = 2;
v_start = 2; a0 = 0.05; eta = 5; % Initial hunger; no neighbours
d_v = 2; d_tau = 6;
hunger_threshold = 0.5; toxin_threshold = 0.1; epsilon = 0.01;
R_search = 4; sigma_weight = 2;
x_domain = -10:0.1:10;
y_domain = -10:0.1:10;
[X_mesh, Y_mesh] = meshgrid(x_domain, y_domain);
centers = [-3*sqrt(3), -3; 0, 6; 3*sqrt(3), -3];
s_peaks = [1, 2, 3]; tau_peaks = [0.8, 1.6, 2.4];
radius_peak = 3; decay = 0.5;
initial_position = [0,0];
domain_limits = [-10,10];
H_threshold = d_v*hunger_threshold;
T_threshold = d_tau*toxin_threshold;

%% Figure settings (defined once)
font_name = 'Times New Roman';
font_size = 18;
food_color = [1,0,0];
avoid_color = [0,0,1];
diffusion_color = [0.45,0.49,0.54];
noise_color = diffusion_color;
actual_color = [0.12,0.12,0.12];
mode_colors = [diffusion_color;food_color;avoid_color]; % actual_code + 1
patch_colors = [0.27,0.46,0.70;0.89,0.47,0.20;0.17,0.63,0.17];
fig1_size = [40,60,1000,600];
fig2_size = [60,90,1100,650];
ax_path_pos = [0.09,0.3,0.38,0.6];
right_x = 0.585; right_width = 0.38; panel_height = 0.175;
right_bottoms = [0.72,0.46,0.20];
strip_bottom = 0.135; strip_height = 0.025;
vector_bottom = 0.32; vector_width = 0.23; vector_height = 0.40;


%% Original food/toxin fields, including the nonzero patch-edge offset
s = zeros(size(X_mesh)); tau = s;
for k = 1:3
    radius = hypot(X_mesh-centers(k,1), Y_mesh-centers(k,2));
    inside = radius <= radius_peak;
    profile = (1 + exp(-decay*radius(inside))) / 2;
    s(inside) = s(inside) + s_peaks(k)*profile;
    tau(inside) = tau(inside) + tau_peaks(k)*profile;
end

%% Fixed centre start and reproducible random increments
input_file = fullfile(repo_dir, 'data', 'agent', 'single_agent_movement_input.csv');
if use_saved_random_input && isfile(input_file)
    input_data = readtable(input_file);
    required = {'step','initial_x','initial_y','z_x','z_y'};
    if ~all(ismember(required, input_data.Properties.VariableNames)) || ...
            height(input_data) ~= num_steps-1
        error('Input CSV needs %d rows: step,initial_x,initial_y,z_x,z_y.',num_steps-1);
    end
    if ~isequal(input_data.step(:), (1:num_steps-1)')
        error('Input CSV steps must be consecutive integers from 1 to %d.',num_steps-1);
    end
    Z = [input_data.z_x, input_data.z_y];
    repeated_xy = [input_data.initial_x, input_data.initial_y];
    if any(~isfinite(repeated_xy(:))) || any(~isfinite(Z(:)))
        error('All input coordinates and random increments must be finite.');
    end
    if any(any(abs(repeated_xy-initial_position) > 1e-12))
        error('The input CSV must use the requested centre start (0,0) on every row.');
    end
    input_source = 'Companion CSV';
else
    rng(random_seed, 'twister');
    Z = randn(2, num_steps-1)';
    input_source = sprintf('MATLAB seed %d', random_seed);
end
if any(initial_position < domain_limits(1)) || any(initial_position > domain_limits(2))
    error('The initial position must be inside [-10,10]^2.');
end

%% Storage: states and the decisions acting between consecutive states
time = (0:num_steps-1)'*dt;
n_decisions = num_steps-1;
positions = zeros(num_steps,2); hunger = zeros(num_steps,1);
food = zeros(num_steps,1); toxin = zeros(num_steps,1); region = zeros(num_steps,1);
H = zeros(n_decisions,1); T = H; Delta = H;
food_gate = false(n_decisions,1); toxin_gate = food_gate;
mode = repmat({'none'},n_decisions,1);        % Eligible behavioural rule
actual_mode = repmat({'diffusion'},n_decisions,1);
actual_code = zeros(n_decisions,1);          % 0 diffusion, 1 food, 2 avoidance
has_directed_movement = false(n_decisions,1);
targets = zeros(n_decisions,2); directions = targets;
drift = targets; noise = sqrt(2*D*dt)*Z;
pre_reflection_displacement = targets; actual_displacement = targets;
reflected = false(n_decisions,1);

positions(1,:) = initial_position; hunger(1) = v_start;
food(1) = interpolate_field(initial_position(1),initial_position(2),s,x_domain,y_domain);
toxin(1) = interpolate_field(initial_position(1),initial_position(2),tau,x_domain,y_domain);
region(1) = patch_region(initial_position,centers,radius_peak);

for n = 1:n_decisions
    old_position = positions(n,:);
    % Original explicit Euler update; u_local=0 for the isolated agent.
    hunger(n+1) = max(0,hunger(n)+dt*(a0-eta*food(n)*hunger(n)));
    H(n) = d_v*hunger(n+1);
    T(n) = d_tau*toxin(n);
    Delta(n) = H(n)-T(n);
    food_gate(n) = H(n) > H_threshold;
    toxin_gate(n) = T(n) > T_threshold;

    target = old_position;
    if food_gate(n) && Delta(n) > epsilon
        mode{n} = 'food';
        target = optimal_target(old_position,s,x_domain,y_domain, ...
            X_mesh,Y_mesh,R_search,sigma_weight,'max');
    elseif toxin_gate(n) && -Delta(n) > epsilon
        mode{n} = 'avoid';
        target = optimal_target(old_position,tau,x_domain,y_domain, ...
            X_mesh,Y_mesh,R_search,sigma_weight,'min');
    end

    targets(n,:) = target;
    vector = target-old_position;
    target_distance = norm(vector);
    if target_distance > 1e-6
        has_directed_movement(n) = true;
        directions(n,:) = vector/target_distance;
        drift(n,:) = p_s*dt*directions(n,:);
        actual_mode{n} = mode{n};
        if strcmp(mode{n},'food'), actual_code(n)=1; else, actual_code(n)=2; end
    end
    pre_reflection_displacement(n,:) = noise(n,:)+drift(n,:);
    proposed_position = old_position+pre_reflection_displacement(n,:);
    next_position = proposed_position;
    % Same sequential reflection rules as the original model.
    if next_position(1)<domain_limits(1), next_position(1)=2*domain_limits(1)-next_position(1); end
    if next_position(1)>domain_limits(2), next_position(1)=2*domain_limits(2)-next_position(1); end
    if next_position(2)<domain_limits(1), next_position(2)=2*domain_limits(1)-next_position(2); end
    if next_position(2)>domain_limits(2), next_position(2)=2*domain_limits(2)-next_position(2); end
    reflected(n) = any(abs(next_position-proposed_position)>1e-12);
    positions(n+1,:) = next_position;
    actual_displacement(n,:) = next_position-old_position;
    food(n+1) = interpolate_field(next_position(1),next_position(2),s,x_domain,y_domain);
    toxin(n+1) = interpolate_field(next_position(1),next_position(2),tau,x_domain,y_domain);
    region(n+1) = patch_region(next_position,centers,radius_peak);
end

%% Representative decisions retained only for the supplementary vectors/data
% Choose the floor midpoint of the longest run of each actual movement class.
% If a class is absent, fill from 25%,50%,75% of the intervals; sort by time.
% The main trajectory and time panels are linked by movement colours only.
event_indices = select_events(actual_code);
event_number = zeros(n_decisions,1);
event_number(event_indices) = (1:numel(event_indices))';

%% Save simulation states and every decision, before plotting
state_table = table(time,positions(:,1),positions(:,2),hunger,food,toxin,region, ...
    'VariableNames',{'time','x','y','hunger','food','toxin','region'});
decision_table = table((1:n_decisions)',time(1:end-1),time(2:end), ...
    positions(1:end-1,1),positions(1:end-1,2),hunger(1:end-1),hunger(2:end), ...
    food(1:end-1),toxin(1:end-1),H,T,Delta,food_gate,toxin_gate, ...
    mode,has_directed_movement,actual_mode,targets(:,1),targets(:,2), ...
    directions(:,1),directions(:,2),drift(:,1),drift(:,2),noise(:,1),noise(:,2), ...
    pre_reflection_displacement(:,1),pre_reflection_displacement(:,2), ...
    actual_displacement(:,1),actual_displacement(:,2), ...
    positions(2:end,1),positions(2:end,2),reflected,event_number, ...
    'VariableNames',{'step','time_start','time_end','x_start','y_start', ...
    'hunger_before','hunger_after','food_local','toxin_local','H','T','Delta', ...
    'food_gate','toxin_gate','mode','has_directed_movement','actual_mode', ...
    'target_x','target_y','direction_x','direction_y','drift_dx','drift_dy', ...
    'noise_dx','noise_dy','pre_reflection_dx','pre_reflection_dy', ...
    'actual_dx','actual_dy','x_end','y_end','reflected','event_number'});
writetable(state_table,fullfile(results_dir,'single_agent_cognition_states.csv'));
writetable(decision_table,fullfile(results_dir,'single_agent_cognition_decisions.csv'));
save(fullfile(results_dir,'single_agent_cognition_results.mat'), ...
    'state_table','decision_table','time','positions','hunger','food','toxin','region', ...
    'H','T','Delta','mode','actual_mode','actual_code','has_directed_movement', ...
    'targets','directions','drift','noise','event_indices', ...
    's','tau','x_domain','y_domain','centers','s_peaks','tau_peaks','radius_peak', ...
    'dt','D','p_s','R_search','sigma_weight','v_start','a0','eta','d_v','d_tau', ...
    'hunger_threshold','toxin_threshold','epsilon','initial_position','Z','input_source');

%% Figure 1: trajectory and decision terms
fig = figure('Name','Movement linked to hunger and toxin', ...
    'NumberTitle','off','Color','w','Units','pixels', ...
    'WindowStyle','normal','Position',fig1_size);
force_light_figure(fig);
ax_path = axes('Parent',fig,'Units','normalized', ...
    'PositionConstraint','innerposition','Position',ax_path_pos);
hold(ax_path,'on');
theta = linspace(0,2*pi,361);
for k = 1:3
    patch(ax_path,centers(k,1)+radius_peak*cos(theta), ...
        centers(k,2)+radius_peak*sin(theta),[0.96,0.96,0.96], ...
        'EdgeColor',[0.68,0.68,0.68],'LineWidth',1.1);
end
% Each actual displacement is coloured by the decision governing that step.
for n = 1:n_decisions
    plot(ax_path,positions(n:n+1,1),positions(n:n+1,2), ...
        'Color',mode_colors(actual_code(n)+1,:),'LineWidth',1);
end
% Mark food-patch centres.
for k = 1:3
    plot(ax_path,centers(k,1),centers(k,2),'*','Color',patch_colors(k,:), ...
        'MarkerSize',22,'LineWidth',1.5,'HandleVisibility','off');
end
h_start = plot(ax_path,positions(1,1),positions(1,2),'o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7,'LineWidth',0.5);
h_end = plot(ax_path,positions(end,1),positions(end,2),'s', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7);
axis(ax_path,'equal'); xlim(ax_path,domain_limits); ylim(ax_path,domain_limits);
xticks(ax_path,domain_limits(1):5:domain_limits(2));
yticks(ax_path,domain_limits(1):5:domain_limits(2));
style_axes(ax_path,font_name,font_size);
h_food = plot(ax_path,nan,nan,'Color',food_color,'LineWidth',2.5);
h_avoid = plot(ax_path,nan,nan,'Color',avoid_color,'LineWidth',2.5);
h_diffusion = plot(ax_path,nan,nan,'Color',diffusion_color,'LineWidth',2.5);
legend(ax_path,[h_food,h_avoid,h_diffusion,h_start,h_end], ...
    {'Food-seeking','Toxin-avoidance','Diffusion only','Start','End'}, ...
    'Location','southoutside','NumColumns',2,'Box','off', ...
    'FontName',font_name,'FontSize',font_size,'TextColor',[0.1,0.1,0.1],'Color','w');

% Values are held on the decision interval [t_n,t_{n+1}). The last step ends
% at 24.95; there is no invented 500th decision at the terminal state.
ax_H = axes('Parent',fig,'Position',[right_x,right_bottoms(1),right_width,panel_height]);
ax_T = axes('Parent',fig,'Position',[right_x,right_bottoms(2),right_width,panel_height]);
ax_delta = axes('Parent',fig,'Position',[right_x,right_bottoms(3),right_width,panel_height]);
upper_HT = max(1,ceil(1.06*max([H;T])));
delta_extent = max(1,ceil(1.08*max(abs(Delta))));
axes_list = [ax_H,ax_T,ax_delta];
values = {H,T,Delta};
titles = {'Weighted hunger: $H_n = d_v v_{n+1}$', ...
          'Weighted local toxin: $T_n = d_\tau \tau(x_n,y_n) $', ...
          'Decision difference: $\Delta_n = H_n - T_n$'};
for j = 1:3
    ax = axes_list(j); hold(ax,'on');
    plot_colored_steps(ax,time,values{j},actual_code,mode_colors);
    xlim(ax,[time(1),time(end)]);
    if j<3, ylim(ax,[0,upper_HT]); else, ylim(ax,[-delta_extent,delta_extent]); end
    h = title(ax,titles{j}, ... 
        'Interpreter','latex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','normal');

    h.Units = 'normalized';
    h.Position(1) = 0;
    h.HorizontalAlignment = 'left';
    style_axes(ax,font_name,font_size);
end
plot(ax_H,[time(1),time(end)],H_threshold*[1,1],'--', ...
    'Color',[0.20,0.20,0.20],'LineWidth',1.1);
plot(ax_T,[time(1),time(end)],T_threshold*[1,1],'--', ...
    'Color',[0.20,0.20,0.20],'LineWidth',1.1);
plot(ax_delta,[time(1),time(end)],epsilon*[1,1],'--', ...
    'Color',[0.20,0.20,0.20],'LineWidth',1.0);
plot(ax_delta,[time(1),time(end)],-epsilon*[1,1],'--', ...
    'Color',[0.20,0.20,0.20],'LineWidth',1.0);
ylabel(ax_H,'$H_n$','Interpreter','latex','FontName',font_name,'FontSize',font_size);
ylabel(ax_T,'$T_n$','Interpreter','latex','FontName',font_name,'FontSize',font_size);
ylabel(ax_delta,'$\Delta_n$','Interpreter','latex','FontName',font_name,'FontSize',font_size);
xlabel(ax_delta,'Decision time, $t_n$','Interpreter','latex','FontName',font_name,'FontSize',font_size);
set(ax_H,'XTickLabel',[]); set(ax_T,'XTickLabel',[]);
linkaxes(axes_list,'x');

%% Action strip below the three decision panels

% Move the shared time-axis labels onto the action strip
strip_ticks = get(ax_delta,'XTick');
set(ax_delta,'XTickLabel',[]);
xlabel(ax_delta,'');

% Align the strip with the right-hand panels
p_delta = get(ax_delta,'Position');

ax_action = axes('Parent',fig,'Units','normalized', ...
    'Position',[p_delta(1),strip_bottom,p_delta(3),strip_height], ...
    'PositionConstraint','innerposition');
hold(ax_action,'on');

% Find consecutive intervals with the same movement mode
run_first = [1; find(diff(actual_code(:)) ~= 0) + 1];
run_last  = [run_first(2:end)-1; numel(actual_code)];

% Draw one coloured rectangle for each interval
for k = 1:numel(run_first)
    n1 = run_first(k);
    n2 = run_last(k);

    t1 = time(n1);
    t2 = time(n2+1);

    patch(ax_action,[t1,t2,t2,t1],[0,0,1,1], ...
        mode_colors(actual_code(n1)+1,:), ...
        'EdgeColor','none', ...
        'HandleVisibility','off');
end

style_axes(ax_action,font_name,font_size);
grid(ax_action,'off');

set(ax_action,'YTick',[],'XTick',strip_ticks);
xlim(ax_action,[time(1),time(end)]);
ylim(ax_action,[0,1]);

xlabel(ax_action,'time, $t_n$', ...
    'Interpreter','latex', ...
    'FontName',font_name,'FontSize',font_size);

linkaxes([axes_list,ax_action],'x');

drawnow;
export_figure(fig,fullfile(figures_dir,'single_agent_movement_link'));

%% Supplement: resolve the directed, random and actual steps at the events
% Tip-to-tail addition: directed step + random step = proposed displacement.
% The actual displacement also includes reflection when a boundary is hit.
% All three panels use the same physical scale; arrows are NOT magnified.
fig_vectors = figure('Name','Representative movement vectors', ...
    'NumberTitle','off','Color','w','Units','pixels', ...
    'WindowStyle','normal','Position',fig2_size);
force_light_figure(fig_vectors);
vector_axes = gobjects(1,numel(event_indices));
mode_names = {'Diffusion only','Food-seeking','Toxin avoidance'};
for j = 1:numel(event_indices)
    n = event_indices(j);
    ax = axes('Parent',fig_vectors, ...
        'Units','normalized', ...
        'Position',[0.095+(j-1)*0.31,vector_bottom,vector_width,vector_height]);
    vector_axes(j) = ax; hold(ax,'on');
    plot(ax,[-0.3,0.3],[0,0],'-','Color',[0.78,0.78,0.78],'LineWidth',0.8);
    plot(ax,[0,0],[-0.3,0.3],'-','Color',[0.78,0.78,0.78],'LineWidth',0.8);
    % Draw actual first so a dashed random arrow remains visible on top
    % when the two coincide in a diffusion-only step.
    quiver(ax,0,0,actual_displacement(n,1),actual_displacement(n,2),0, ...
        'Color',actual_color,'LineWidth',1.7,'MaxHeadSize',0.8);
    quiver(ax,drift(n,1),drift(n,2),noise(n,1),noise(n,2),0, ...
        'Color',noise_color,'LineStyle','--','LineWidth',1.7,'MaxHeadSize',0.8);
    if has_directed_movement(n)
        quiver(ax,0,0,drift(n,1),drift(n,2),0, ...
            'Color',mode_colors(actual_code(n)+1,:),'LineWidth',1.7,'MaxHeadSize',0.8);
    end
    plot(ax,0,0,'o','Color','k','MarkerFaceColor','w','MarkerSize',5,'LineWidth',1.0);
    axis(ax,'equal'); xlim(ax,[-0.3,0.3]); ylim(ax,[-0.3,0.3]);
    xticks(ax,[-0.2,0,0.2]); yticks(ax,[-0.2,0,0.2]);
    
    xlabel(ax,'$x$ displacement','Interpreter','latex','FontName',font_name,'FontSize',font_size);
    if j==1, ylabel(ax,'$y$ displacement','Interpreter','latex','FontName',font_name,'FontSize',font_size); end
    % Hunger threshold comparison
    if H(n) > H_threshold
        H_text = '$H_n > d_v v_h$';
    else
        H_text = '$H_n \leq d_v v_h$';
    end

    % Toxin threshold comparison
    if T(n) > T_threshold
        T_text = '$T_n > d_\tau \tau_h$';
    else
        T_text = '$T_n \leq d_\tau \tau_h$';
    end

    % Decision difference comparison
    if Delta(n) > epsilon
        delta_text = '$\Delta_n > \epsilon$';

    elseif Delta(n) < -epsilon
        delta_text = '$\Delta_n < -\epsilon$';

    else
        delta_text = '$|\Delta_n| \leq \epsilon$';
    end

    title(ax, ...
        {sprintf('%s $\\mid\\; t=%.2f$', ...
        mode_names{actual_code(n)+1}, time(n)), ...
        sprintf('%s, %s, %s', ...
        H_text, T_text, delta_text)}, ...
        'Interpreter','latex', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'FontWeight','normal');

    style_axes(ax,font_name,font_size);
end
linkaxes(vector_axes,'xy');



% One common legend, rather than repeating four entries inside each panel.
%% Figure 2: common legend with arrows

% Two rows, two entries per row
% Normalized coordinates relative to the full figure

legend_x = [0.16, 0.55];
legend_y = [0.115, 0.065];

arrow_len = 0.045;

legend_labels = { ...
    'Directed: food-seeking', ...
    'Random displacement', ...
    'Directed: toxin avoidance', ...
    'Actual displacement'};

legend_colors = [ ...
    food_color;
    noise_color;
    avoid_color;
    actual_color];

legend_styles = {'-','--','-','-'};

for k = 1:4

    col = mod(k-1,2) + 1;
    row = floor((k-1)/2) + 1;

    x0 = legend_x(col);
    y0 = legend_y(row);

    annotation(fig_vectors,'arrow', ...
        [x0, x0+arrow_len], [y0,y0], ...
        'Color',legend_colors(k,:), ...
        'LineStyle',legend_styles{k}, ...
        'LineWidth',1.8, ...
        'HeadLength',4, ...
        'HeadWidth',7);

    annotation(fig_vectors,'textbox', ...
        [x0+arrow_len+0.01,y0-0.022,0.28,0.045], ...
        'String',legend_labels{k}, ...
        'Interpreter','none', ...
        'FontName',font_name, ...
        'FontSize',font_size, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', ...
        'EdgeColor','none', ...
        'FitBoxToText','off');
end

drawnow;
export_figure(fig_vectors,fullfile(figures_dir,'single_agent_movement_vectors'));

fprintf('Saved %d states and %d decisions. Input: %s\n',num_steps,n_decisions,input_source);
fprintf('Actual steps: food=%d; avoidance=%d; diffusion only=%d.\n', ...
    sum(actual_code==1),sum(actual_code==2),sum(actual_code==0));
fprintf('Eligible rules without directed movement: %d.\n', ...
    sum(~strcmp(mode,'none') & ~has_directed_movement));
disp(decision_table(event_indices,{'step','time_start','H','T','Delta','mode','actual_mode'}));

%% Local helpers: source-faithful interpolation and weighted target selection
function val = interpolate_field(x,y,field,x_grid,y_grid)
ix = max(1,min(length(x_grid)-1,find(x_grid<=x,1,'last')));
iy = max(1,min(length(y_grid)-1,find(y_grid<=y,1,'last')));
if isempty(ix) || isempty(iy), val=0; return; end
wx = (x-x_grid(ix))/(x_grid(ix+1)-x_grid(ix));
wy = (y-y_grid(iy))/(y_grid(iy+1)-y_grid(iy));
val = (1-wx)*(1-wy)*field(iy,ix) + wx*(1-wy)*field(iy,ix+1) ...
    + (1-wx)*wy*field(iy+1,ix) + wx*wy*field(iy+1,ix+1);
end

function target = optimal_target(position,field,x_grid,y_grid,X_mesh,Y_mesh,R,sigma,optimize_type)
[~,ix] = min(abs(x_grid-position(1))); [~,iy] = min(abs(y_grid-position(2)));
rx = ceil(R/(x_grid(2)-x_grid(1))); ry = ceil(R/(y_grid(2)-y_grid(1)));
columns = max(1,ix-rx):min(length(x_grid),ix+rx);
rows = max(1,iy-ry):min(length(y_grid),iy+ry);
subX = X_mesh(rows,columns); subY = Y_mesh(rows,columns);
subfield = field(rows,columns);
distance2 = (subX-position(1)).^2+(subY-position(2)).^2;
within = distance2<=R^2;
target = position;
if ~any(within(:)), return; end
candidateX = subX(within); candidateY = subY(within); candidatefield = subfield(within);
% The source checks raw field uniformity BEFORE applying perception weights.
if max(candidatefield)-min(candidatefield)<1e-10, return; end
weight = exp(-distance2(within)/(2*max(sigma,1e-12)^2));
if strcmp(optimize_type,'max')
    objective = weight.*candidatefield; optimum = max(objective);
else
    objective = candidatefield./max(weight,1e-12); optimum = min(objective);
end
candidates = find(abs(objective-optimum)<1e-12);
if numel(candidates)>1
    tied_distance2 = (candidateX(candidates)-position(1)).^2 ...
        +(candidateY(candidates)-position(2)).^2;
    [~,nearest] = min(tied_distance2); chosen = candidates(nearest);
else
    chosen = candidates(1);
end
target = [candidateX(chosen),candidateY(chosen)];
end

function region = patch_region(position,centers,radius_peak)
[distance,region] = min(hypot(position(1)-centers(:,1),position(2)-centers(:,2)));
if distance>radius_peak, region=0; end
end

function indices = select_events(code)
indices = [];
for class = [1,2,0]
    mask = code==class;
    starts = find(diff([false;mask;false])==1);
    ends = find(diff([false;mask;false])==-1)-1;
    if ~isempty(starts)
        [~,longest] = max(ends-starts+1);
        indices(end+1) = floor((starts(longest)+ends(longest))/2); %#ok<AGROW>
    end
end
fill_indices = floor((numel(code)-1)*[0.25,0.5,0.75])+1;
for index = fill_indices
    if numel(indices)>=3, break; end
    if ~ismember(index,indices), indices(end+1)=index; end %#ok<AGROW>
end
indices = sort(indices(:));
end

function plot_colored_steps(ax,time,values,code,colors)
% A gray stair underlay connects vertical changes between different modes.
stairs(ax,time,[values;values(end)],'Color',[0.65,0.65,0.65],'LineWidth',0.8);
starts = [1;find(diff(code)~=0)+1]; ends = [starts(2:end)-1;numel(code)];
for k = 1:numel(starts)
    intervals = starts(k):ends(k);
    stairs(ax,[time(intervals);time(ends(k)+1)], ...
        [values(intervals);values(ends(k))], ...
        'Color',colors(code(starts(k))+1,:),'LineWidth',1.8);
end
end

function style_axes(ax,font_name,font_size)
set(ax,'FontName',font_name,'FontSize',font_size,'LineWidth',1.1, ...
    'Box','on','TickDir','out','Layer','top','Color','w', ...
    'XTickLabelRotation',0,'TickLabelInterpreter','latex', ...
    'XColor',[0.1,0.1,0.1],'YColor',[0.1,0.1,0.1], ...
    'ZColor',[0.1,0.1,0.1],'GridColor',[0.5,0.5,0.5]);
set([ax.XLabel,ax.YLabel,ax.ZLabel,ax.Title],'Color',[0.1,0.1,0.1]);
grid(ax,'on'); ax.GridAlpha=0.18;
end

function force_light_figure(fig)
% Keep this figure readable under dark desktop themes or custom defaults.
% Theme is available in R2025a+; explicit colours also support older releases.
if isprop(fig,'Theme')
    fig.Theme = 'light';
end
set(fig,'Color','w');
for ax = reshape(findall(fig,'Type','axes'),1,[])
    set(ax,'Color','w','XColor',[0.1,0.1,0.1], ...
        'YColor',[0.1,0.1,0.1],'ZColor',[0.1,0.1,0.1], ...
        'GridColor',[0.5,0.5,0.5]);
end
% Labels, titles and annotations need their own text colours.
set(findall(fig,'Type','text'),'Color',[0.1,0.1,0.1]);
set(findall(fig,'Type','textboxshape'),'Color',[0.1,0.1,0.1]);
set(findall(fig,'Type','legend'),'TextColor',[0.1,0.1,0.1],'Color','w');
end

function export_figure(fig,basename)
force_light_figure(fig);
drawnow;
set(fig,'PaperPositionMode','auto');
if exist('exportgraphics','file')==2
    exportgraphics(fig,[basename,'.eps'],'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig,[basename,'.png'],'Resolution',300,'BackgroundColor','white');
    exportgraphics(fig,[basename,'.pdf'],'ContentType','vector','BackgroundColor','white');
else
    print(fig,[basename,'.eps'],'-depsc','-painters');
    print(fig,[basename,'.png'],'-dpng','-r300');
    print(fig,[basename,'.pdf'],'-dpdf','-painters','-bestfit');
end
end

