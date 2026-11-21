% =========================================================================
% 2D AGENT-BASED MODEL: FOOD-TOXIN TRADEOFF WITH HUNGER AND COGNITION
%
% Movement based on MAXIMUM perceived food w_R*s / MINIMUM perceived toxin tau/w_R
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

tic;

%% GLOBAL PARAMETERS
num_steps = 500;       % can reduce for faster testing (e.g., 200)
dt = 0.05;
num_agents = 500;      % can reduce to 100 for testing

% Domain
x_min = -10; x_max = 10; dx = 0.1;
y_min = -10; y_max = 10; dy = 0.1;
x_domain = x_min:dx:x_max;
y_domain = y_min:dy:y_max;
[X_mesh, Y_mesh] = meshgrid(x_domain, y_domain);

% Hunger parameters
v_start = 10;  % Default initial hunger; rerun for 0, 2, and 10 to compare.
v_start_override = getenv('AGENT_INITIAL_HUNGER');
if ~isempty(v_start_override)
    v_start = str2double(v_start_override);
    if ~isfinite(v_start) || v_start < 0
        error('AGENT_INITIAL_HUNGER must be a finite, nonnegative number.');
    end
end
a0 = 0.05;
eta = 5;               % food-mediated hunger reduction rate
beta = 5;              % crowding strength in 1 + beta*u_local

% Hunger-toxin decision parameters
hunger_threshold = 0.5; % v_h
toxin_threshold  = 0.1; % tau_h
epsilon = 0.01;          % decision margin from Table: epsilon

% Perception parameters (matched to the continuum formulation)
R_l = 0.5;              % local perception radius
R_n = 4.0;              % nonlocal perception radius
sigma_factor = 0.5;          % fixed Gaussian spatial scale

% Movement
D = 2;                 % diffusion
p_s = 2;               % step size for directed movement

% Death/reproduction (not fully used)
r = 0.3;
d1 = 0.5;
d_v = 2;
d_tau = 6;

% Experimental conditions
conditions = {'Random', 'NoCognition', 'Cognition'};

%% CREATE FOOD (s) AND TOXIN (tau) FIELDS (3 Gaussian peaks)
centers = {[-3*sqrt(3), -3], [0, 6], [3*sqrt(3), -3]};
s_peaks = [1, 2, 3];
tau_peaks = [0.8, 1.6, 2.4];
toxin_mode = getenv('AGENT_TOXIN_MODE');
if isempty(toxin_mode), toxin_mode = 'Toxin'; end
if strcmpi(toxin_mode, 'noToxin')
    toxin_mode = 'noToxin';
    tau_peaks = [0, 0, 0];
elseif strcmpi(toxin_mode, 'Toxin')
    toxin_mode = 'Toxin';
else
    error('AGENT_TOXIN_MODE must be Toxin or noToxin.');
end
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

%% RUN MULTIPLE REPETITIONS FOR EACH CONDITION
num_reps = 50;   % number of simulations per condition
fprintf('Running %d repetitions per condition...\n', num_reps);

% Preallocate: conditions x regions x repetitions
all_agent_counts = zeros(length(conditions), 3, num_reps);

for cond_idx = 1:length(conditions)
    condition = conditions{cond_idx};
    fprintf('\n========================================\n');
    fprintf('Condition: %s\n', condition);
    fprintf('========================================\n');

    for rep = 1:num_reps
        % Run one simulation
        results = run_simulation(condition, ...
            num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
            s, tau, v_start, a0, eta, beta, d_v, d_tau, ...
            hunger_threshold, toxin_threshold, epsilon, D, p_s, r, d1, ...
            R_l, R_n, sigma_factor, centers, radius_peak);

        % Count agents in each region at final time step
        radius = radius_peak;  % 3
        region_counts = zeros(1, 3);
        for i = 1:num_agents
            if results.alive(i, end)
                final_x = results.positions(i, 1, end);
                final_y = results.positions(i, 2, end);
                % Find closest region centre
                dists = zeros(1,3);
                for k = 1:3
                    dists(k) = norm([final_x - centers{k}(1), final_y - centers{k}(2)]);
                end
                [min_dist, closest] = min(dists);
                if min_dist <= radius
                    region_counts(closest) = region_counts(closest) + 1;
                end
            end
        end
        all_agent_counts(cond_idx, :, rep) = region_counts;

        if mod(rep, 10) == 0
            fprintf('  Rep %d/%d done\n', rep, num_reps);
        end
    end
end

%% COMPUTE STATISTICS ACROSS REPETITIONS
mean_counts = mean(all_agent_counts, 3);
min_counts = min(all_agent_counts, [], 3);
max_counts = max(all_agent_counts, [], 3);



% Name the file after its actual parameter: the original script used v_start=10
% but wrote 'v0_...', which could silently mislabel the results.
output_name = sprintf('v%g_%s_agent_counts_%dreps.mat', ...
    v_start, toxin_mode, num_reps);
save(fullfile(results_dir, output_name), 'all_agent_counts', 'conditions', ...
    'num_agents', 'v_start', 'num_reps', 'toxin_mode');
fprintf('Saved %s\n', fullfile(results_dir, output_name));

%% Print statistics to command window
fprintf('\n=== Summary across %d repetitions ===\n', num_reps);
for cond_idx = 1:length(conditions)
    fprintf('\n%s:\n', conditions{cond_idx});
    for reg = 1:3
        fprintf('  Region %d: mean = %.1f, min = %d, max = %d\n', ...
            reg, mean_counts(cond_idx,reg), min_counts(cond_idx,reg), max_counts(cond_idx,reg));
    end
end

toc;

%% ==================== HELPER FUNCTIONS ====================
% (all the function definitions – interpolate_field, compute_local_density,
% find_optimal_point_weighted, run_simulation – go here, exactly as you wrote them.)


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
    x_grid, y_grid, X_mesh, Y_mesh, R, sigma, optimize_type, weight_mode)
% Finds an optimal perceived location within radius R.
%
% Food seeking:
%   perceived food = w_R(z;x) * s(z), then maximize.
%
% Toxin avoidance:
%   perceived toxin = tau(z) / w_R(z;x), then minimize.
%
% The Gaussian weight is
%   w_R(z;x) = exp(-|z-x|^2/(2*sigma^2)).
%
% If the environmental field itself is uniform in the search region,
% return the current position so there is no artificial directed motion.

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
    best_x = x;
    best_y = y;
    return;
end

% Gaussian perception weight
weights = exp(-dist2 / (2*sigma^2));

% Keep only points inside the circular perception region
sub_X = sub_X(within);
sub_Y = sub_Y(within);
sub_field = sub_field(within);
weights = weights(within);

% Uniform environmental cue -> no directed movement
if max(sub_field) - min(sub_field) < 1e-10
    best_x = x;
    best_y = y;
    return;
end

% Construct perceived objective
switch weight_mode
    case 'multiply'
        % Food: maximize w_R * s
        perceived_field = sub_field .* weights;
    case 'divide'
        % Toxin: minimize tau / w_R
        perceived_field = sub_field ./ max(weights, 1e-12);
    otherwise
        error('weight_mode must be ''multiply'' or ''divide''');
end

% Optimize perceived objective
switch optimize_type
    case 'max'
        opt_val = max(perceived_field);
        candidates = find(abs(perceived_field - opt_val) <= 1e-12);
    case 'min'
        opt_val = min(perceived_field);
        candidates = find(abs(perceived_field - opt_val) <= 1e-12);
    otherwise
        error('optimize_type must be ''max'' or ''min''');
end

% If several locations tie, choose the nearest one (largest Gaussian weight)
[~, nearest_idx] = max(weights(candidates));
idx_opt = candidates(nearest_idx);

best_x = sub_X(idx_opt);
best_y = sub_Y(idx_opt);
end

%% ==================== SIMULATION FUNCTION ====================

function results = run_simulation(condition, ...
    num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
    s, tau, v_start, a0, eta, beta, d_v, d_tau, ...
    hunger_threshold, toxin_threshold, epsilon, D, p_s, r, d1, ...
    R_l, R_n, sigma_factor, centers, radius_peak)

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
        R_search = R_l;
        sigma_weight = sigma_factor * R_search;
        current_d_v = d_v;

    case 'Cognition'
        use_max_movement = true;
        R_search = R_n;
        sigma_weight = sigma_factor * R_search;
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
        % dv_i/dt = a0 - eta*s(X_i)/(1 + beta*u_i)
        dv_dt = a0 - eta * s_val / (1 + beta * u_local);
        v_new = v + dt * dv_dt;
        v_new = max(0, v_new);
        hunger_levels(i,t) = v_new;

        % --- Movement ---
        % 1. Diffusion
        displacement = sqrt(2*D*dt) * randn(2,1);

        % 2. Directed movement (only if condition uses it)
        if use_max_movement
            % Competing internal/environmental pressures
            hunger_pressure = current_d_v * v_new;
            toxin_pressure  = d_tau * tau_val;
            pressure_diff   = hunger_pressure - toxin_pressure;

            if v_new > hunger_threshold && pressure_diff >= epsilon
                % Hunger dominates toxin by at least epsilon:
                % z_s^R in argmax w_R(z;x)*s(z).
                % Food target maximizes w_R(z;x) * s(z).
                [best_x, best_y] = find_optimal_point_weighted(x, y, s, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, ...
                    'max', 'multiply');

            elseif tau_val > toxin_threshold && pressure_diff <= -epsilon
                % Toxin dominates hunger by at least epsilon:
                % z_tau^R in argmin tau(z)/w_R(z;x).
                % Toxin target minimizes tau(z) / w_R(z;x).
                [best_x, best_y] = find_optimal_point_weighted(x, y, tau, ...
                    x_domain, y_domain, X_mesh, Y_mesh, R_search, sigma_weight, ...
                    'min', 'divide');

            else
                % |hunger_pressure - toxin_pressure| < epsilon,
                % or neither threshold is active -> diffusion only.
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

