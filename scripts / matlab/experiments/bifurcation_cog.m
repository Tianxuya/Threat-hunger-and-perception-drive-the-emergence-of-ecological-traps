%% BIFURCATION CALCULATION
% R = 0   : Random movement
% R > 0   : Cognition with nonlocal perception radius R_n = R
%
% This version is aligned with the parameter definitions and helper
% functions in the current 2D agent-based food-toxin model.

clear; close all; clc;

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

%% ==================== BIFURCATION PARAMETERS ====================

cognition_radii = 0:0.5:10;
num_reps = 50;

%% ==================== GLOBAL PARAMETERS ====================

num_steps = 500;
dt = 0.05;
num_agents = 500;

% Domain
x_min = -10; x_max = 10; dx = 0.1;
y_min = -10; y_max = 10; dy = 0.1;

x_domain = x_min:dx:x_max;
y_domain = y_min:dy:y_max;
[X_mesh, Y_mesh] = meshgrid(x_domain, y_domain);

% Hunger parameters
v_start = 0;
v_start_override = getenv('AGENT_INITIAL_HUNGER');
if ~isempty(v_start_override)
    v_start = str2double(v_start_override);
    if ~isfinite(v_start) || v_start < 0
        error('AGENT_INITIAL_HUNGER must be a finite, nonnegative number.');
    end
end
a0 = 0.05;
eta = 5;
beta = 5;

% Hunger-toxin decision parameters
hunger_threshold = 0.5;
toxin_threshold  = 0.1;
epsilon = 0.01;

% Perception parameters
R_l = 0.5;              % local perception radius
sigma_factor = 0.5;     % sigma = sigma_factor * R_search

% Movement
D = 2;
p_s = 2;

% Death/reproduction parameters (currently not actively used)
r = 0.3;
d1 = 0.5;
d_v = 2;
d_tau = 6;

%% ==================== FOOD AND TOXIN FIELDS ====================

centers = {[-3*sqrt(3), -3], ...
           [0, 6], ...
           [3*sqrt(3), -3]};

s_peaks   = [1, 2, 3];
tau_peaks = [0.8, 1.6, 2.4];

radius_peak = 3;
decay = 0.5;

s = zeros(size(X_mesh));
tau = zeros(size(X_mesh));

for i = 1:size(X_mesh,1)
    for j = 1:size(X_mesh,2)
        for k = 1:3

            dist = norm([X_mesh(i,j) - centers{k}(1), ...
                         Y_mesh(i,j) - centers{k}(2)]);

            if dist <= radius_peak

                exp_factor = exp(-decay * dist);

                s(i,j) = s(i,j) + ...
                    s_peaks(k)/2 + s_peaks(k)*exp_factor/2;

                tau(i,j) = tau(i,j) + ...
                    tau_peaks(k)/2 + tau_peaks(k)*exp_factor/2;
            end
        end
    end
end

s = max(0, s);
tau = max(0, tau);

%% ==================== PREALLOCATE BIFURCATION RESULTS ====================

num_radii = length(cognition_radii);

% repetitions x regions x cognition-radius
all_counts = zeros(num_reps, 3, num_radii);

mean_counts = zeros(num_radii, 3);
std_counts  = zeros(num_radii, 3);
min_counts  = zeros(num_radii, 3);
max_counts  = zeros(num_radii, 3);

%% ==================== BIFURCATION LOOP ====================

for r_idx = 1:num_radii

    R = cognition_radii(r_idx);

    fprintf('\n========================================\n');
    fprintf('Perception radius R = %.1f\n', R);
    fprintf('========================================\n');

    % R = 0 is the random baseline.
    % Every R > 0 uses the Cognition condition, with R_n = R.
    if R == 0
        condition = 'Random';

        % R_n is unused by Random because directed movement is disabled.
        R_n = R_l;
    else
        condition = 'Cognition';
        R_n = R;
    end

    counts_this_radius = zeros(num_reps, 3);

    for rep = 1:num_reps

        results = run_simulation(condition, ...
            num_steps, dt, num_agents, ...
            x_domain, y_domain, X_mesh, Y_mesh, ...
            s, tau, ...
            v_start, a0, eta, beta, ...
            d_v, d_tau, ...
            hunger_threshold, toxin_threshold, epsilon, ...
            D, p_s, r, d1, ...
            R_l, R_n, sigma_factor, ...
            centers, radius_peak);

        %% Count final agents in the three resource patches

        region_counts = zeros(1,3);

        for i = 1:num_agents

            if results.alive(i,end)

                final_x = results.positions(i,1,end);
                final_y = results.positions(i,2,end);

                % Assign an agent to the closest patch centre,
                % provided that it lies within radius_peak.
                dists = zeros(1,3);

                for k = 1:3
                    dists(k) = norm([final_x - centers{k}(1), ...
                                     final_y - centers{k}(2)]);
                end

                [min_dist, closest_region] = min(dists);

                if min_dist <= radius_peak
                    region_counts(closest_region) = ...
                        region_counts(closest_region) + 1;
                end
            end
        end

        counts_this_radius(rep,:) = region_counts;
        all_counts(rep,:,r_idx) = region_counts;

        if mod(rep,10) == 0
            fprintf('  Rep %d/%d done\n', rep, num_reps);
        end
    end

    %% Statistics across repetitions

    mean_counts(r_idx,:) = mean(counts_this_radius, 1);
    std_counts(r_idx,:)  = std(counts_this_radius, 0, 1);
    min_counts(r_idx,:)  = min(counts_this_radius, [], 1);
    max_counts(r_idx,:)  = max(counts_this_radius, [], 1);

    fprintf(['  Low food:    %.1f +/- %.1f\n' ...
             '  Medium food: %.1f +/- %.1f\n' ...
             '  High food:   %.1f +/- %.1f\n'], ...
             mean_counts(r_idx,1), std_counts(r_idx,1), ...
             mean_counts(r_idx,2), std_counts(r_idx,2), ...
             mean_counts(r_idx,3), std_counts(r_idx,3));
end

%% ==================== SAVE RESULTS ====================

output_name = sprintf('bifurcation_all_patches_v%g.mat', v_start);
save(fullfile(results_dir, output_name), ...
    'cognition_radii', ...
    'all_counts', ...
    'mean_counts', ...
    'std_counts', ...
    'min_counts', ...
    'max_counts', ...
    'num_reps', ...
    'num_steps', ...
    'dt', ...
    'num_agents', ...
    'v_start', ...
    'a0', ...
    'eta', ...
    'beta', ...
    'hunger_threshold', ...
    'toxin_threshold', ...
    'epsilon', ...
    'R_l', ...
    'sigma_factor', ...
    'D', ...
    'p_s', ...
    'r', ...
    'd1', ...
    'd_v', ...
    'd_tau', ...
    'centers', ...
    's_peaks', ...
    'tau_peaks', ...
    'radius_peak', ...
    'decay');

fprintf('\nData saved to %s\n', fullfile(results_dir, output_name));

%% ==================== PLOT ====================

figure;

subplot(3,1,1);
fill_x = [cognition_radii, fliplr(cognition_radii)];
fill_y = [mean_counts(:,1)' + std_counts(:,1)', ...
          fliplr(mean_counts(:,1)' - std_counts(:,1)')];
fill(fill_x, fill_y, [0.8 0.8 0.8], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(cognition_radii, mean_counts(:,1), 'LineWidth', 2);
ylabel('Agents');
title('Low-food patch');
grid on;

subplot(3,1,2);
fill_y = [mean_counts(:,2)' + std_counts(:,2)', ...
          fliplr(mean_counts(:,2)' - std_counts(:,2)')];
fill(fill_x, fill_y, [0.8 0.8 0.8], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(cognition_radii, mean_counts(:,2), 'LineWidth', 2);
ylabel('Agents');
title('Medium-food patch');
grid on;

subplot(3,1,3);
fill_y = [mean_counts(:,3)' + std_counts(:,3)', ...
          fliplr(mean_counts(:,3)' - std_counts(:,3)')];
fill(fill_x, fill_y, [0.8 0.8 0.8], ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(cognition_radii, mean_counts(:,3), 'LineWidth', 2);
xlabel('Cognition radius R');
ylabel('Agents');
title('High-food patch');
grid on;

exportgraphics(gcf, fullfile(figures_dir, sprintf('bifurcation_all_patches_v%g.png', v_start)), 'Resolution', 300);

toc;

%% ==================== HELPER FUNCTIONS ====================

function val = interpolate_field(x, y, field, x_grid, y_grid)
% Bilinear interpolation with boundary clamping

nx = length(x_grid);
ny = length(y_grid);

ix = max(1, min(nx-1, find(x_grid <= x, 1, 'last')));
iy = max(1, min(ny-1, find(y_grid <= y, 1, 'last')));

if isempty(ix) || isempty(iy) || ix == nx || iy == ny
    val = 0;
    return;
end

x1 = x_grid(ix);
x2 = x_grid(ix+1);
y1 = y_grid(iy);
y2 = y_grid(iy+1);

wx = (x - x1) / (x2 - x1);
wy = (y - y1) / (y2 - y1);

f11 = field(iy, ix);
f21 = field(iy, ix+1);
f12 = field(iy+1, ix);
f22 = field(iy+1, ix+1);

val = (1-wx)*(1-wy)*f11 + ...
      wx*(1-wy)*f21 + ...
      (1-wx)*wy*f12 + ...
      wx*wy*f22;
end


function u_local = compute_local_density(positions, x, y, radius)
% Count agents within the local density radius.

if isempty(positions)
    u_local = 0;
    return;
end

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
% Gaussian perception weight:
%   w_R(z;x) = exp(-|z-x|^2/(2*sigma^2)).
%
% If the environmental field itself is uniform in the perception region,
% return the current position to avoid artificial directed movement.

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

% Prevent division by zero if a very small R is ever tested.
sigma = max(sigma, 1e-12);

% Gaussian perception weight
weights = exp(-dist2 / (2*sigma^2));

% Keep only candidates within the circular perception region
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

% If several points tie, choose the nearest one.
% Because Gaussian weight decreases with distance, this is equivalent
% to choosing the tied candidate with largest weight.
[~, nearest_idx] = max(weights(candidates));
idx_opt = candidates(nearest_idx);

best_x = sub_X(idx_opt);
best_y = sub_Y(idx_opt);
end


function results = run_simulation(condition, ...
    num_steps, dt, num_agents, x_domain, y_domain, X_mesh, Y_mesh, ...
    s, tau, v_start, a0, eta, beta, d_v, d_tau, ...
    hunger_threshold, toxin_threshold, epsilon, D, p_s, r, d1, ...
    R_l, R_n, sigma_factor, centers, radius_peak)

fprintf('Running condition: %s\n', condition);

%% Set perception/movement mode

switch condition

    case 'Random'
        use_max_movement = false;
        R_search = 0;
        sigma_weight = 0;
        current_d_v = 0;

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

%% Initialise agents on radius-1 circle around origin

positions = zeros(num_agents, 2, num_steps);
hunger_levels = zeros(num_agents, num_steps);
alive = true(num_agents, num_steps);

theta = linspace(0, 2*pi, num_agents+1)';

positions(:,1,1) = cos(theta(1:end-1));
positions(:,2,1) = sin(theta(1:end-1));

hunger_levels(:,1) = v_start;

%% Statistics

region_counts = zeros(3, num_steps);
mean_hunger = zeros(1, num_steps);
total_agents = zeros(1, num_steps);

total_agents(1) = num_agents;
mean_hunger(1) = v_start;

%% Time loop

for t = 2:num_steps

    current_pos = positions(:,:,t-1);
    current_hunger = hunger_levels(:,t-1);
    current_alive = alive(:,t-1);

    alive_count = sum(current_alive);
    total_agents(t) = alive_count;

    if alive_count == 0
        break;
    end

    for i = 1:num_agents

        if ~current_alive(i)
            continue;
        end

        x = current_pos(i,1);
        y = current_pos(i,2);
        v = current_hunger(i);

        %% Environmental values at current position

        s_val = interpolate_field(x, y, s, x_domain, y_domain);
        tau_val = interpolate_field(x, y, tau, x_domain, y_domain);

        %% Local density

        u_local = compute_local_density( ...
            current_pos(current_alive,:), x, y, 1.0);

        %% Hunger dynamics
        % dv_i/dt = a0 - eta*s(X_i)/(1 + beta*u_i)

        dv_dt = a0 - eta * s_val / (1 + beta * u_local);

        v_new = v + dt * dv_dt;
        v_new = max(0, v_new);

        hunger_levels(i,t) = v_new;

        %% Movement: diffusion

        displacement = sqrt(2*D*dt) * randn(2,1);

        %% Directed movement

        if use_max_movement

            hunger_pressure = current_d_v * v_new;
            toxin_pressure  = d_tau * tau_val;
            pressure_diff   = hunger_pressure - toxin_pressure;

            if v_new > hunger_threshold && pressure_diff >= epsilon

                % Hunger dominates toxin:
                % maximize w_R(z;x) * s(z)

                [best_x, best_y] = find_optimal_point_weighted( ...
                    x, y, s, ...
                    x_domain, y_domain, X_mesh, Y_mesh, ...
                    R_search, sigma_weight, ...
                    'max', 'multiply');

            elseif tau_val > toxin_threshold && pressure_diff <= -epsilon

                % Toxin dominates hunger:
                % minimize tau(z) / w_R(z;x)

                [best_x, best_y] = find_optimal_point_weighted( ...
                    x, y, tau, ...
                    x_domain, y_domain, X_mesh, Y_mesh, ...
                    R_search, sigma_weight, ...
                    'min', 'divide');

            else

                % Neither pressure dominates by epsilon:
                % diffusion only.

                best_x = x;
                best_y = y;
            end

            dir_vec = [best_x - x; best_y - y];
            dist_to_target = norm(dir_vec);

            if dist_to_target > 1e-6

                dir_unit = dir_vec / dist_to_target;

                displacement = displacement + ...
                    p_s * dir_unit * dt;
            end
        end

        %% Position update with reflecting boundaries

        x_new = x + displacement(1);
        y_new = y + displacement(2);

        if x_new < x_domain(1)
            x_new = 2*x_domain(1) - x_new;
        end

        if x_new > x_domain(end)
            x_new = 2*x_domain(end) - x_new;
        end

        if y_new < y_domain(1)
            y_new = 2*y_domain(1) - y_new;
        end

        if y_new > y_domain(end)
            y_new = 2*y_domain(end) - y_new;
        end

        positions(i,:,t) = [x_new, y_new];

        %% Region assignment for statistics

        dist_region1 = norm([x_new - centers{1}(1), ...
                             y_new - centers{1}(2)]);

        dist_region2 = norm([x_new - centers{2}(1), ...
                             y_new - centers{2}(2)]);

        dist_region3 = norm([x_new - centers{3}(1), ...
                             y_new - centers{3}(2)]);

        [min_region_dist, region] = ...
            min([dist_region1, dist_region2, dist_region3]);

        if min_region_dist <= radius_peak
            region_counts(region,t) = region_counts(region,t) + 1;
        end

        %% Reproduction/death
        % Currently disabled, matching the present model.

        alive(i,t) = true;
    end

    % Previously dead agents remain dead
    alive(:,t) = alive(:,t) & current_alive;

    if sum(current_alive) > 0
        mean_hunger(t) = mean(hunger_levels(current_alive,t));
    else
        mean_hunger(t) = 0;
    end
end

%% Package results

results = struct();

results.condition = condition;
results.positions = positions;
results.hunger_levels = hunger_levels;
results.alive = alive;
results.region_counts = region_counts;
results.mean_hunger = mean_hunger;
results.total_agents = total_agents;

end
