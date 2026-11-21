% =========================================================================
% 2D AGENT-BASED MODEL: FOOD-TOXIN TRADEOFF WITH HUNGER AND COGNITION
%
% Movement based on MAXIMUM WEIGHTED food / MINIMUM WEIGHTED toxin
% within a search radius, using a Gaussian kernel centered on the agent.
% =========================================================================

clear all; close all; clc;

% Resolve outputs from this file, independent of the current directory.
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(fileparts(script_dir)));
results_dir = getenv('AGENT_RESULTS_DIR');
if isempty(results_dir), results_dir = fullfile(repo_dir, 'results', 'agent'); end
figures_dir = getenv('AGENT_FIGURES_DIR');
if isempty(figures_dir), figures_dir = fullfile(repo_dir, 'figures', 'agent'); end
if ~isfolder(results_dir), mkdir(results_dir); end
if ~isfolder(figures_dir), mkdir(figures_dir); end


%% GLOBAL PARAMETERS
num_steps = 500;       % Reduced for testing (increase if needed)
dt = 0.05;
num_agents = 500;

% Domain
x_min = -10; x_max = 10; dx = 0.1;
y_min = -10; y_max = 10; dy = 0.1;
x_domain = x_min:dx:x_max;
y_domain = y_min:dy:y_max;
[X_mesh, Y_mesh] = meshgrid(x_domain, y_domain);

% Hunger parameters
v_start = 10; 
a0 = 0.05;
b = 0.2;

% Cognition parameters
cognition_radius = 4.0;
R_nocog = 0.5;         % Tiny radius for NoCognition
sigma_factor = 0.5;    % weight sigma = R * sigma_factor

% Movement
D = 2;                 % diffusion
p_s = 2;               % step size for directed movement

% Death/reproduction
r = 0.8;
d1 = 0.03;
d_v = 2;
d_tau = 6;

% Experimental conditions
conditions = {'Random', 'NoCognition', 'Cognition'};

%% CREATE FOOD (s) AND TOXIN (tau) FIELDS (3 Gaussian peaks)
% Three centres with different food‑to‑toxin ratios
centers = {[-3*sqrt(3), -3], [0, 6], [3*sqrt(3), -3]};
s_peaks = [1, 2, 3];
 % tau_peaks = [0., 0, 0];
 tau_peaks = [0.8, 1.6, 2.4];
radius_peak = 3;
decay = 0.5;

s = zeros(size(X_mesh));
tau = zeros(size(X_mesh));

for i = 1:size(X_mesh,1)
    for j = 1:size(X_mesh,2)
        for k = 1:3
            dist = norm([X_mesh(i,j)-centers{k}(1), Y_mesh(i,j)-centers{k}(2)]);
            if dist <= radius_peak
                exp_factor = exp(-decay * dist);
                s(i,j) = s(i,j) + s_peaks(k)/2 + s_peaks(k)*exp_factor/2;
                tau(i,j) = tau(i,j) + tau_peaks(k)/2 + tau_peaks(k)*exp_factor/2;
            end
        end
    end
end
s = max(0, s);
tau = max(0, tau);

%% RUN ALL CONDITIONS
all_results = struct();
for cond_idx = 1:length(conditions)
    condition = conditions{cond_idx};
    results = run_simulation(condition, ...
        num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
        s, tau, v_start, a0, b, d_v, d_tau, D, p_s, r, d1, ...
        cognition_radius, R_nocog, sigma_factor, centers, radius_peak);
    all_results.(condition) = results;
end

% Save the simulated fields and trajectories for stand-alone plotting scripts.
save(fullfile(results_dir, 'food_toxin_max_determined_reaction.mat'), ...
    'all_results', 'conditions', 'centers', 'num_agents', 'radius_peak', ...
    'x_domain', 'y_domain', 'X_mesh', 'Y_mesh', 's', 'tau');

%% ==================== HELPER FUNCTIONS ====================

function val = interpolate_field(x, y, field, x_grid, y_grid)
% Bilinear interpolation with boundary clamping
nx = length(x_grid); ny = length(y_grid);
ix = max(1, min(nx-1, find(x_grid <= x, 1, 'last')));
iy = max(1, min(ny-1, find(y_grid <= y, 1, 'last')));
if isempty(ix) || isempty(iy) || ix == nx || iy == ny
    val = 0; return;
end
x1 = x_grid(ix); x2 = x_grid(ix+1);
y1 = y_grid(iy); y2 = y_grid(iy+1);
wx = (x - x1) / (x2 - x1);
wy = (y - y1) / (y2 - y1);
f11 = field(iy, ix);   f21 = field(iy, ix+1);
f12 = field(iy+1, ix); f22 = field(iy+1, ix+1);
val = (1-wx)*(1-wy)*f11 + wx*(1-wy)*f21 + (1-wx)*wy*f12 + wx*wy*f22;
end

function u_local = compute_local_density(positions, x, y, radius)
% Count neighbours within radius (excluding self)
if isempty(positions); u_local = 0; return; end
dx = positions(:,1) - x;
dy = positions(:,2) - y;
u_local = sum(dx.^2 + dy.^2 < radius^2) - 1;   % subtract self
u_local = max(u_local, 0);
end

function [best_x, best_y] = find_optimal_point_weighted(x, y, field, ...
    x_grid, y_grid, X_mesh, Y_mesh, R, sigma, optimize_type)
% Finds point within radius R that maximizes (or minimizes) the
% Gaussian‑weighted field. If all weighted values are equal (within 1e-10),
% returns (x,y) to signal "no directed movement".

[~, idx_x] = min(abs(x_grid - x));
[~, idx_y] = min(abs(y_grid - y));
dx_step = x_grid(2)-x_grid(1);
dy_step = y_grid(2)-y_grid(1);
rad_steps_x = ceil(R / dx_step);
rad_steps_y = ceil(R / dy_step);
i_min = max(1, idx_x - rad_steps_x);
i_max = min(length(x_grid), idx_x + rad_steps_x);
j_min = max(1, idx_y - rad_steps_y);
j_max = min(length(y_grid), idx_y + rad_steps_y);

sub_X = X_mesh(j_min:j_max, i_min:i_max);
sub_Y = Y_mesh(j_min:j_max, i_min:i_max);
sub_field = field(j_min:j_max, i_min:i_max);

dist2 = (sub_X - x).^2 + (sub_Y - y).^2;
within = dist2 <= R^2;

if ~any(within(:))
    best_x = x; best_y = y;
    return;
end

% Apply Gaussian weight
weights = exp(-dist2 / (2*sigma^2));
weighted_field = sub_field .* weights;

% Keep only points inside the circle
sub_X = sub_X(within);
sub_Y = sub_Y(within);
weighted_field = weighted_field(within);

% Check if all weighted values are equal (uniform field)
max_val = max(weighted_field);
min_val = min(weighted_field);
if abs(max_val - min_val) < 1e-10
    % Uniform field → no directed movement
    best_x = x;
    best_y = y;
    return;
end

switch optimize_type
    case 'max'
        [~, idx_opt] = max(weighted_field);
    case 'min'
        [~, idx_opt] = min(weighted_field);
    otherwise
        error('optimize_type must be ''max'' or ''min''');
end
best_x = sub_X(idx_opt);
best_y = sub_Y(idx_opt);
end

%% ==================== SIMULATION FUNCTION ====================

function results = run_simulation(condition, ...
    num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
    s, tau, v_start, a0, b, d_v, d_tau, D, p_s, r, d1, ...
    cognition_radius, R_nocog, sigma_factor, centers, radius_peak)

fprintf('Running condition: %s\n', condition);

% Set parameters based on condition
switch condition
    case 'Random'
        use_max_movement = false;
        R_search = 0;
        sigma_weight = 0;
        current_d_v = 0;
    case 'NoCognition'
        use_max_movement = true;
        R_search = R_nocog;
        sigma_weight = R_nocog * sigma_factor;
        current_d_v = d_v;
    case 'Cognition'
        use_max_movement = true;
        R_search = cognition_radius;
        sigma_weight = cognition_radius * sigma_factor;
        current_d_v = d_v;
    otherwise
        error('Unknown condition');
end

max_agents = 2000;   % maximum capacity (avoid unbounded growth)
positions = zeros(max_agents, 2);
hunger = zeros(max_agents, 1);
alive = false(max_agents, 1);

% Initialise agents on a circle of radius 1 around origin
theta = linspace(0, 2*pi, num_agents+1)';
init_x = cos(theta(1:end-1));
init_y = sin(theta(1:end-1));
for i = 1:num_agents
    positions(i,1) = init_x(i);
    positions(i,2) = init_y(i);
    hunger(i) = v_start;
    alive(i) = true;
end
agent_count = num_agents;

% Storage for statistics (time series)
region_counts = zeros(3, num_steps);
mean_hunger = zeros(1, num_steps);
total_agents = zeros(1, num_steps);

epsilon = 0.1;   % small constant to avoid division by zero in birth rate

for t = 2:num_steps
    % Compute current positions of active agents
    current_pos = positions(alive, :);
    current_hunger = hunger(alive);
    alive_indices = find(alive);
    N = length(alive_indices);
    total_agents(t) = N;
    if N == 0
        break;
    end
    
    % Record mean hunger
    mean_hunger(t) = mean(current_hunger);
    
    for idx = 1:N
        i = alive_indices(idx);
        x = positions(i,1);
        y = positions(i,2);
        v = hunger(i);
        
        % Interpolate environment at current position
        s_val = interpolate_field(x, y, s, x_domain, y_domain);
        tau_val = interpolate_field(x, y, tau, x_domain, y_domain);
        
        % Local density (number of neighbours within radius 1)
        u_local = compute_local_density(current_pos, x, y, 1.0);
        
        % --- Hunger dynamics ---
        dv_dt = a0 - s_val / (b + u_local);
        v_new = v + dt * dv_dt;
        v_new = max(0, v_new);
        hunger(i) = v_new;
        
        % --- Movement ---
        displacement = sqrt(2*D*dt) * randn(2,1);
        
        if use_max_movement
            hunger_threshold = 0.5;
            toxin_threshold = 0.1;
            if v_new > hunger_threshold && current_d_v * v_new > d_tau * tau_val
                [best_x, best_y] = find_optimal_point_weighted(x, y, s, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, 'max');
            elseif tau_val > toxin_threshold && current_d_v * v_new <= d_tau * tau_val
                [best_x, best_y] = find_optimal_point_weighted(x, y, tau, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, 'min');
            else
                best_x = x; best_y = y;
            end
            dir_vec = [best_x - x; best_y - y];
            dist = norm(dir_vec);
            if dist > 1e-6
                dir_unit = dir_vec / dist;
                displacement = displacement + p_s * dir_unit * dt;
            end
        end
        
        x_new = x + displacement(1);
        y_new = y + displacement(2);
        % Reflection on boundaries
        if x_new < x_domain(1); x_new = 2*x_domain(1) - x_new; end
        if x_new > x_domain(end); x_new = 2*x_domain(end) - x_new; end
        if y_new < y_domain(1); y_new = 2*y_domain(1) - y_new; end
        if y_new > y_domain(end); y_new = 2*y_domain(end) - y_new; end
        positions(i,1) = x_new;
        positions(i,2) = y_new;
        
        % --- Birth and Death (independent probabilities) ---
        % Birth rate: logistic term depending on local density and food
        birth_rate = r;
        % Death rate: baseline * hunger factor * toxin factor
        death_rate = d1 * (1 + d_v * v_new) * (1 + d_tau * tau_val);
        
        p_birth = min(birth_rate * dt, 1);
        p_death = min(death_rate * dt, 1);
        
        % Check death first (if the agent dies, it cannot reproduce)
        if rand() < p_death
            alive(i) = false;
        elseif rand() < p_birth
            % Create a new agent at the same location (clone)
            if agent_count < max_agents
                agent_count = agent_count + 1;
                positions(agent_count,:) = [x_new, y_new];
                hunger(agent_count) = v_new;   % offspring inherits hunger
                alive(agent_count) = true;
            end
        end
    end
    
    % Update region counts for this time step (based on all alive agents)
    cnt = zeros(1,3);
    for i = find(alive)'
        x = positions(i,1); y = positions(i,2);
        for k = 1:3
            dist = norm([x - centers{k}(1), y - centers{k}(2)]);
            if dist <= radius_peak
                cnt(k) = cnt(k) + 1;
                break;
            end
        end
    end
    region_counts(:,t) = cnt';
    
    if mod(t,200) == 0
        fprintf('  Step %d, Agents: %d, Mean hunger: %.3f\n', t, sum(alive), mean_hunger(t));
    end
end

% Package results
results = struct();
results.condition = condition;
results.positions = positions(1:agent_count, :);
results.hunger_levels = hunger(1:agent_count);
results.alive = alive(1:agent_count);
results.region_counts = region_counts;
results.mean_hunger = mean_hunger;
results.total_agents = total_agents;
fprintf('Condition %s completed. Final agents: %d\n\n', condition, sum(alive));
end