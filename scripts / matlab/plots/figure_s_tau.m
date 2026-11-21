clear;
paths = project_paths();
clc;
close all;

%% ============================================================
%  SPATIAL DOMAIN
% ============================================================
outdir = paths.agent_figures;
if ~isfolder(outdir), mkdir(outdir); end



x = linspace(-10, 10, 301);
y = linspace(-10, 10, 301);

[X, Y] = meshgrid(x, y);


%% ============================================================
%  PATCH CENTERS
% ============================================================

center1 = [-3*sqrt(3), -3];
center2 = [0, 6];
center3 = [3*sqrt(3), -3];


%% ============================================================
%  FOOD AND TOXIN PEAK VALUES
% ============================================================

s1_peak = 1;
s2_peak = 2;
s3_peak = 3;

tau1_peak = 0.8;
tau2_peak = 1.6;
tau3_peak = 2.4;


%% ============================================================
%  PATCH PARAMETERS
% ============================================================

radius = 3;
decay_rate = 0.5;

s = zeros(size(X));
tau = zeros(size(X));


%% ============================================================
%  CONSTRUCT FOOD AND TOXIN DISTRIBUTIONS
%
%  For patch k:
%
%  s_k(x,y)
%  =
%  (s_k^max / 2)
%  [1 + exp(-alpha r_k)]
%  for r_k <= R
%
%  tau_k(x,y)
%  =
%  (tau_k^max / 2)
%  [1 + exp(-alpha r_k)]
%  for r_k <= R
%
% ============================================================

for i = 1:size(X,1)

    for j = 1:size(X,2)

        %% Distance to each patch center

        dist1 = norm( ...
            [X(i,j)-center1(1), ...
             Y(i,j)-center1(2)] );

        dist2 = norm( ...
            [X(i,j)-center2(1), ...
             Y(i,j)-center2(2)] );

        dist3 = norm( ...
            [X(i,j)-center3(1), ...
             Y(i,j)-center3(2)] );


        %% ---------------- PATCH 1 ----------------

        if dist1 <= radius

            exp_factor1 = exp(-decay_rate * dist1);

            s(i,j) = s(i,j) ...
                + s1_peak/2 ...
                + s1_peak*exp_factor1/2;

            tau(i,j) = tau(i,j) ...
                + tau1_peak/2 ...
                + tau1_peak*exp_factor1/2;

        end


        %% ---------------- PATCH 2 ----------------

        if dist2 <= radius

            exp_factor2 = exp(-decay_rate * dist2);

            s(i,j) = s(i,j) ...
                + s2_peak/2 ...
                + s2_peak*exp_factor2/2;

            tau(i,j) = tau(i,j) ...
                + tau2_peak/2 ...
                + tau2_peak*exp_factor2/2;

        end


        %% ---------------- PATCH 3 ----------------

        if dist3 <= radius

            exp_factor3 = exp(-decay_rate * dist3);

            s(i,j) = s(i,j) ...
                + s3_peak/2 ...
                + s3_peak*exp_factor3/2;

            tau(i,j) = tau(i,j) ...
                + tau3_peak/2 ...
                + tau3_peak*exp_factor3/2;

        end

    end

end


%% Ensure non-negative values

s = max(0, s);
tau = max(0, tau);


%% ============================================================
%  CUSTOM COLORMAPS
% ============================================================

nColors = 256;

%% Food: very light green -> dark green

food_cmap = [
    linspace(0.95, 0.08, nColors)', ...
    linspace(1.00, 0.50, nColors)', ...
    linspace(0.95, 0.15, nColors)'
];


%% Toxin: very light blue -> dark blue

toxin_cmap = [
    linspace(0.95, 0.05, nColors)', ...
    linspace(0.98, 0.40, nColors)', ...
    linspace(1.00, 0.80, nColors)'
];


%% ============================================================
%  FIGURE 1: FOOD AND TOXIN DISTRIBUTIONS
% ============================================================

fig=figure( ...
    'Color', 'w', ...
    'Position', [5, 5, 800, 400],Visible='off');


%% ============================================================
%  FOOD DISTRIBUTION
% ============================================================

subplot(1,2,1);

contourf( ...
    X, Y, s, ...
    20, ...
    'LineColor', 'none');

colormap(gca, food_cmap);

cb1 = colorbar;
cb1.FontSize = 14;
cb1.Label.String = 'Food availability';
cb1.Label.FontSize = 15;

hold on;


%% Add contour lines

contour( ...
    X, Y, s, ...
    10, ...
    'k', ...
    'LineWidth', 0.5);


%% Add L/M/H labels

text( ...
    center1(1), center1(2), ...
    'L', ...
    'Color', 'k', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);

text( ...
    center2(1), center2(2), ...
    'M', ...
    'Color', 'k', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);

text( ...
    center3(1), center3(2), ...
    'H', ...
    'Color', 'w', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);


hold off;


title( ...
    'Food Distribution $s(x,y)$', ...
    'Interpreter', 'latex', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

xlabel( ...
    '$x$', ...
    'Interpreter', 'latex', ...
    'FontSize', 20);

ylabel( ...
    '$y$', ...
    'Interpreter', 'latex', ...
    'FontSize', 20);

axis equal;

xlim([-10 10]);
ylim([-10 10]);

grid on;

set(gca, ...
    'FontSize', 20, ...
    'LineWidth', 1.0, ...
    'TickDir', 'out', ...
    'GridAlpha', 0.15);


%% ============================================================
%  TOXIN DISTRIBUTION
% ============================================================

subplot(1,2,2);

contourf( ...
    X, Y, tau, ...
    20, ...
    'LineColor', 'none');

colormap(gca, toxin_cmap);

cb2 = colorbar;
cb2.FontSize = 14;
cb2.Label.String = 'Toxin level';
cb2.Label.FontSize = 15;

hold on;


%% Add contour lines

contour( ...
    X, Y, tau, ...
    10, ...
    'k', ...
    'LineWidth', 0.5);


%% Add L/M/H labels

text( ...
    center1(1), center1(2), ...
    'L', ...
    'Color', 'k', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);

text( ...
    center2(1), center2(2), ...
    'M', ...
    'Color', 'k', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);

text( ...
    center3(1), center3(2), ...
    'H', ...
    'Color', 'w', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 17);


hold off;


title( ...
    'Toxin Distribution $\tau(x,y)$', ...
    'Interpreter', 'latex', ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

xlabel( ...
    '$x$', ...
    'Interpreter', 'latex', ...
    'FontSize', 20);

ylabel( ...
    '$y$', ...
    'Interpreter', 'latex', ...
    'FontSize', 20);

axis equal;

xlim([-10 10]);
ylim([-10 10]);

grid on;

set(gca, ...
    'FontSize', 20, ...
    'LineWidth', 1.0, ...
    'TickDir', 'out', ...
    'GridAlpha', 0.15);


exportgraphics(fig, ...
    fullfile(outdir, "food and toxin distribution.eps"), ...
    "ContentType", "vector");
exportgraphics(fig, ...
    fullfile(outdir, "food and toxin distribution.jpg"));
