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
num_agents = 1; %200;

% Domain
x_min = -10; x_max = 10; dx = 0.1;
y_min = -10; y_max = 10; dy = 0.1;
x_domain = x_min:dx:x_max;
y_domain = y_min:dy:y_max;
[X_mesh, Y_mesh] = meshgrid(x_domain, y_domain);

% Hunger parameters
v_start = 10;
a0 = 0.05;
eta = 5;               % equivalent to 1/b when b = 0.2
beta = 5;              % equivalent to 1/b when b = 0.2

% Cognition parameters
cognition_radius = 4.0;
R_nocog = 0.5;         % Tiny radius for NoCognition
sigma_factor = cognition_radius * 0.5;    % weight sigma = R * sigma_factor

% Movement
D = 1;                 % diffusion
p_s = 2;               % step size for directed movement

% Death/reproduction (not fully used, but keep for completeness)
r = 0.3;
d1 = 0.03;
d_v = 2;
d_tau = 60;

% Decision margin
epsilon = 0.01;

% Experimental conditions
conditions = {'Random', 'NoCognition', 'Cognition'};

%% CREATE FOOD (s) AND TOXIN (tau) FIELDS (3 Gaussian peaks)
% Three centres with different food‑to‑toxin ratios
centers = {[-3*sqrt(3), -3], [0, 6], [3*sqrt(3), -3]};
s_peaks = [1, 2, 3];
% tau_peaks = [0 0 0];
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
        s, tau, v_start, a0, eta, beta, d_v, d_tau, epsilon, D, p_s, r, d1, ...
        cognition_radius, R_nocog, sigma_factor, centers, radius_peak);
    all_results.(condition) = results;
end

% Save the simulated fields and trajectories for stand-alone plotting scripts.
save(fullfile(results_dir, 'food_toxin_max_determined.mat'), ...
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
% Count neighbours within radius (excluding self if needed, but fine)
if isempty(positions); u_local = 0; return; end
dx = positions(:,1) - x;
dy = positions(:,2) - y;
u_local = sum(dx.^2 + dy.^2 < radius^2);
end

function [best_x, best_y] = find_optimal_point_weighted(x, y, field, ...
    x_grid, y_grid, X_mesh, Y_mesh, R, sigma, optimize_type)

% Find the preferred point within perception radius R.
%
% Food seeking ('max'):
%       maximize  w_R(z;x) * s(z)
%
% Toxin avoidance ('min'):
%       minimize  tau(z) / w_R(z;x)
%
% where
%       w_R(z;x) = exp(-|z-x|^2/(2*sigma^2)).
%
% If the objective is effectively uniform, the current position is returned,
% corresponding to no directed movement.

% Find grid point nearest to current agent position
[~, idx_x] = min(abs(x_grid - x));
[~, idx_y] = min(abs(y_grid - y));

dx_step = x_grid(2) - x_grid(1);
dy_step = y_grid(2) - y_grid(1);

rad_steps_x = ceil(R / dx_step);
rad_steps_y = ceil(R / dy_step);

i_min = max(1, idx_x - rad_steps_x);
i_max = min(length(x_grid), idx_x + rad_steps_x);

j_min = max(1, idx_y - rad_steps_y);
j_max = min(length(y_grid), idx_y + rad_steps_y);

% Extract local search region
sub_X = X_mesh(j_min:j_max, i_min:i_max);
sub_Y = Y_mesh(j_min:j_max, i_min:i_max);
sub_field = field(j_min:j_max, i_min:i_max);

% Distance from current position
dist2 = (sub_X - x).^2 + (sub_Y - y).^2;

% Restrict candidates to perception circle
within = dist2 <= R^2;

if ~any(within(:))
    best_x = x;
    best_y = y;
    return;
end

% Gaussian perception weight
sigma = max(sigma, 1e-12);

weights = exp(-dist2 / (2*sigma^2));

% Keep only locations within perception radius
candidate_X = sub_X(within);
candidate_Y = sub_Y(within);
candidate_field = sub_field(within);
candidate_weights = weights(within);

switch optimize_type

    case 'max'
        % Food-seeking score:
        % Phi_s = w_R * s
        objective = candidate_field .* candidate_weights;

    case 'min'
        % Toxin-avoidance score:
        % Phi_tau = tau / w_R
        %
        % Avoid division by extremely small Gaussian weights.
        objective = candidate_field ./ max(candidate_weights, 1e-12);

    otherwise
        error('optimize_type must be ''max'' or ''min''');
end

% If all candidate objective values are effectively identical,
% there is no preferred direction.
if max(objective) - min(objective) < 1e-10
    best_x = x;
    best_y = y;
    return;
end

% Determine optimum
switch optimize_type

    case 'max'
        optimum_value = max(objective);
        tied = abs(objective - optimum_value) < 1e-12;

    case 'min'
        optimum_value = min(objective);
        tied = abs(objective - optimum_value) < 1e-12;
end

% If several locations have the same optimum, choose the nearest one.
% This avoids an arbitrary direction caused by MATLAB's first-index rule.
candidate_indices = find(tied);

if length(candidate_indices) > 1

    tied_dist2 = ...
        (candidate_X(candidate_indices) - x).^2 + ...
        (candidate_Y(candidate_indices) - y).^2;

    [~, nearest_idx] = min(tied_dist2);

    idx_opt = candidate_indices(nearest_idx);

else

    idx_opt = candidate_indices(1);

end

best_x = candidate_X(idx_opt);
best_y = candidate_Y(idx_opt);

end

function results = run_simulation(condition, ...
    num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
    s, tau, v_start, a0, eta, beta, d_v, d_tau, epsilon, D, p_s, r, d1, ...
    cognition_radius, R_nocog, sigma_factor, centers, radius_peak)

fprintf('Running condition: %s\n', condition);

% Set parameters based on condition
switch condition
    case 'Random'
        use_max_movement = false;
        R_search = 0;          % not used
        sigma_weight = 0;      % not used
        current_d_v = 0;       % no hunger effect on movement
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

% Initialise agents on a circle of radius 1 around origin
positions = zeros(num_agents, 2, num_steps);
hunger_levels = zeros(num_agents, num_steps);
alive = true(num_agents, num_steps);
theta = linspace(0, 2*pi, num_agents+1)';
positions(:,1,1) = cos(theta(1:end-1));
positions(:,2,1) = sin(theta(1:end-1));
hunger_levels(:,1) = v_start;

% Statistics
region_counts = zeros(3, num_steps);
mean_hunger = zeros(1, num_steps);
total_agents = zeros(1, num_steps);

for t = 2:num_steps
    current_pos = positions(:,:,t-1);
    current_hunger = hunger_levels(:,t-1);
    current_alive = alive(:,t-1);
    alive_count = sum(current_alive);
    total_agents(t) = alive_count;
    if alive_count == 0; break; end

    for i = 1:num_agents
        if ~current_alive(i); continue; end

        x = current_pos(i,1);
        y = current_pos(i,2);
        v = current_hunger(i);

        % Interpolate food and toxin at current position
        s_val = interpolate_field(x, y, s, x_domain, y_domain);
        tau_val = interpolate_field(x, y, tau, x_domain, y_domain);

        % Local density (for hunger update)
        u_local = compute_local_density(current_pos(current_alive,:), x, y, 1.0);

        % Update hunger:
        % dv/dt = a0 - eta*s/(1 + beta*u)
        % With eta = beta = 5, this is equivalent to s/(0.2 + u).
        dv_dt = a0 - eta * s_val / (1 + beta * u_local);
        v_new = v + dt * dv_dt;
        v_new = max(0, v_new);
        hunger_levels(i,t) = v_new;

        % --- Movement ---
        % 1. Diffusion
        displacement = sqrt(2*D*dt) * randn(2,1);

        % 2. Directed movement (only if condition uses it)
        hunger_threshold = 0.5;   % adjust based on your hunger scale
        toxin_threshold  = 0.1;   % adjust based on your toxin scale

        if use_max_movement
            if v_new > hunger_threshold && current_d_v * v_new - d_tau * tau_val > epsilon
                % Hungry enough -> seek food
                [best_x, best_y] = find_optimal_point_weighted(x, y, s, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, 'max');
            elseif tau_val > toxin_threshold && d_tau * tau_val - current_d_v * v_new > epsilon
                % Toxin significant and not hungry enough -> avoid toxin
                [best_x, best_y] = find_optimal_point_weighted(x, y, tau, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, 'min');
            else
                % Neither hunger nor toxin are motivating -> no directed movement
                best_x = x;
                best_y = y;
            end
            dir_vec = [best_x - x; best_y - y];
            dist = norm(dir_vec);
            if dist > 1e-6
                dir_unit = dir_vec / dist;
                displacement = displacement + p_s * dir_unit * dt;
            end
        end

        % Update position with reflection on boundaries
        x_new = x + displacement(1);
        y_new = y + displacement(2);
        if x_new < x_domain(1); x_new = 2*x_domain(1) - x_new; end
        if x_new > x_domain(end); x_new = 2*x_domain(end) - x_new; end
        if y_new < y_domain(1); y_new = 2*y_domain(1) - y_new; end
        if y_new > y_domain(end); y_new = 2*y_domain(end) - y_new; end

        positions(i,:,t) = [x_new, y_new];

        % Region assignment for statistics
        d1 = norm([x_new - centers{1}(1), y_new - centers{1}(2)]);
        d2 = norm([x_new - centers{2}(1), y_new - centers{2}(2)]);
        d3 = norm([x_new - centers{3}(1), y_new - centers{3}(2)]);
        [~, region] = min([d1, d2, d3]);
        if min([d1,d2,d3]) <= radius_peak
            region_counts(region, t) = region_counts(region, t) + 1;
        end

        % Reproduction/death (simplified: always stay alive)
        alive(i,t) = true;
    end

    % Update alive status: if an agent died in previous step, keep dead
    alive(:,t) = alive(:,t) & current_alive;
    if sum(current_alive) > 0
        mean_hunger(t) = mean(hunger_levels(current_alive,t));
    else
        mean_hunger(t) = 0;
    end

    if mod(t, 200) == 0
        fprintf('  Step %d/%d, Alive: %d, Mean hunger: %.3f\n', ...
            t, num_steps, sum(alive(:,t)), mean_hunger(t));
    end
end

% Package results
results = struct();
results.condition = condition;
results.positions = positions;
results.hunger_levels = hunger_levels;
results.alive = alive;
results.region_counts = region_counts;
results.mean_hunger = mean_hunger;
results.total_agents = total_agents;

fprintf('Condition %s completed. Final agents: %d\n\n', ...
    condition, sum(alive(:,end)));
end

