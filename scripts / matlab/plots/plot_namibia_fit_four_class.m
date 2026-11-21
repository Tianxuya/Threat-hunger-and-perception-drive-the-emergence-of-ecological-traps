clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================

% Choose a scenario: 'full', 'nofood', 'notoxin', or 'uniform_low_food'.
scenario = 'full';
indir = fullfile(paths.namibia_results, scenario);
outdir = paths.namibia_figures;


%% ================= READ DATA =================

PDEgrid = readtable(fullfile(indir, "pde_mortality_grid.csv"));
mort    = readtable(fullfile(indir, "mortality_event_comparison.csv"));

lon   = PDEgrid.x;
lat   = PDEgrid.y;
death = PDEgrid.cumulative_total_mortality_mass;


%% ================= CLEAN MODEL DATA =================

% valid = ...
%     isfinite(lon) & ...
%     isfinite(lat) & ...
%     isfinite(death);
% 
% lon   = lon(valid);
% lat   = lat(valid);
% death = death(valid);


%% ================= CLEAN OBSERVED DATA =================

if all(ismember(["X", "Y", "N_Dead"], mort.Properties.VariableNames))

    validObs = ...
        isfinite(mort.X) & ...
        isfinite(mort.Y) & ...
        isfinite(mort.N_Dead);

    mort = mort(validObs, :);

    % Remove observations with 80 deaths or more
    mort = mort(mort.N_Dead < 80, :);

    fprintf("Maximum remaining observed N_Dead = %.0f\n", ...
        max(mort.N_Dead));

end


%% ================= CLIP EXTREME MODEL VALUES =================

% Clip upper 1% of modeled mortality values
vmax = prctile(death, 99);

death_clip = min(death, vmax);


%% ================= INTERPOLATION GRID =================

nx = 300;
ny = 300;

xq = linspace(min(lon), max(lon), nx);
yq = linspace(min(lat), max(lat), ny);

[Xq, Yq] = meshgrid(xq, yq);


%% ================= INTERPOLATE =================

F = scatteredInterpolant( ...
    lon, ...
    lat, ...
    death_clip, ...
    "natural", ...
    "none");

Zq = F(Xq, Yq);


%% ================= PRESERVE NO-DATA REGIONS =================

idx = knnsearch( ...
    [lon lat], ...
    [Xq(:) Yq(:)]);

nearestDist = sqrt( ...
    (Xq(:) - lon(idx)).^2 + ...
    (Yq(:) - lat(idx)).^2);

nearestDist = reshape(nearestDist, size(Xq));

% Smaller value -> preserve larger gaps / holes
distThresh = 0.12;

Zq(nearestDist > distThresh) = NaN;


%% ================= CREATE FOUR MODEL CLASSES =================

validZ = isfinite(Zq);

vals = Zq(validZ);

% Quartile boundaries
q1 = quantile(vals, 0.25);
q2 = quantile(vals, 0.50);
q3 = quantile(vals, 0.75);

fprintf("\nFour-class thresholds:\n");
fprintf("Q1 = %.6g\n", q1);
fprintf("Q2 = %.6g\n", q2);
fprintf("Q3 = %.6g\n\n", q3);


Zclass = nan(size(Zq));

% Class 1: Very low
Zclass(validZ & Zq <= q1) = 1;

% Class 2: Low
Zclass(validZ & Zq > q1 & Zq <= q2) = 2;

% Class 3: Medium
Zclass(validZ & Zq > q2 & Zq <= q3) = 3;

% Class 4: High
Zclass(validZ & Zq > q3) = 4;


%% ================= FIGURE SETTINGS =================

fontSize = 18;
fontName = "Arial";

xlimMap = [11.99489 21.0];
ylimMap = [-29.0 -16.9059];

fig = figure( ...
    "Color", "w", ...
    "Position", [100 100 900 900],"Visible", "off");

ax = axes(fig);


%% ================= PLOT FOUR MODEL CLASSES =================

% contourf(ax, ...
%     Xq, ...
%     Yq, ...
%     Zclass, ...
%     [1 2 3 4 5], ...
%     "LineColor", "none");

hold(ax, "on");


%% ================= FOUR-COLOR COLORMAP =================

cmap = [
    1.00 1.00 0.75   % Very low
    1.00 0.85 0.45   % Low
    1.00 0.60 0.30   % Medium
    0.80 0.10 0.10   % High
];

% colormap(ax, cmap);
% 
% clim(ax, [1 5]);


%% ================= COLORBAR =================

% cb = colorbar(ax);
% 
% cb.Ticks = [1.5 2.5 3.5 4.5];
% 
% cb.TickLabels = { ...
%     "Very low", ...
%     "Low", ...
%     "Medium", ...
%     "High"};
% 
% cb.FontSize = fontSize;
% cb.FontName = fontName;
% cb.TickLength = 0;


%% ================= OBSERVED MORTALITY DATA =================

if all(ismember(["X", "Y", "N_Dead"], mort.Properties.VariableNames))

    % ---------------------------------------------------------
    % Exact observed-death classes:
    %
    % 1--3
    % 4--6
    % 7--10
    % 11--16
    % 17--24
    % 25--33
    % 34--50
    % 51--79
    % ---------------------------------------------------------

    labels = { ...
        "1--3", ...
        "4--6", ...
        "7--10", ...
        "11--16", ...
        "17--24", ...
        "25--33", ...
        "34--50", ...
        ">50"};

    % Marker sizes for the 8 classes
    sizeClass = [ ...
        35 ...
        55 ...
        80 ...
        110 ...
        150 ...
        200 ...
        270 ...
        350];

    % Initialize marker sizes
    obs_size = nan(size(mort.N_Dead));


    % ---------------------------------------------------------
    % Assign exact classes
    % ---------------------------------------------------------

    obs_size(mort.N_Dead >= 1  & mort.N_Dead <= 3)  = sizeClass(1);

    obs_size(mort.N_Dead >= 4  & mort.N_Dead <= 6)  = sizeClass(2);

    obs_size(mort.N_Dead >= 7  & mort.N_Dead <= 10) = sizeClass(3);

    obs_size(mort.N_Dead >= 11 & mort.N_Dead <= 16) = sizeClass(4);

    obs_size(mort.N_Dead >= 17 & mort.N_Dead <= 24) = sizeClass(5);

    obs_size(mort.N_Dead >= 25 & mort.N_Dead <= 33) = sizeClass(6);

    obs_size(mort.N_Dead >= 34 & mort.N_Dead <= 50) = sizeClass(7);

    obs_size(mort.N_Dead >= 51 & mort.N_Dead < 80)  = sizeClass(8);


    % Only plot observations successfully assigned to a class
    validPlot = isfinite(obs_size);


    %% ================= PLOT OBSERVED POINTS =================

    scatter(ax, ...
        mort.X(validPlot), ...
        mort.Y(validPlot), ...
        obs_size(validPlot), ...
        "MarkerEdgeColor", "w", ...
        "MarkerFaceColor", "k", ...
        "MarkerFaceAlpha", 0.75, ...
        "LineWidth", 1.0);


    %% ================= MANUAL BUBBLE LEGEND =================

    xLeg = 12.5;
    yTop = -23.8;
    dy   = 0.58;


    % text(ax, ...
    %     xLeg - 0.15, ...
    %     yTop + 0.55, ...
    %     "Observed deaths", ...
    %     "FontSize", 14, ...
    %     "FontName", fontName, ...
    %     "HorizontalAlignment", "left");


    for k = 1:numel(sizeClass)

        yk = yTop - (k-1)*dy;


        % Bubble
        scatter(ax, ...
            xLeg, ...
            yk, ...
            sizeClass(k), ...
            "MarkerEdgeColor", "w", ...
            "MarkerFaceColor", "k", ...
            "MarkerFaceAlpha", 0.75, ...
            "LineWidth", 1.0);


        % Label
        % text(ax, ...
        %     xLeg + 0.35, ...
        %     yk, ...
        %     labels{k}, ...
        %     "FontSize", 14, ...
        %     "FontName", fontName, ...
        %     "VerticalAlignment", "middle", ...
        %     "HorizontalAlignment", "left");

    end

end


% ================= AXIS LABELS =================

% xlabel(ax, ...
%     "Longitude", ...
%     "FontSize", fontSize, ...
%     "FontName", fontName);
% 
% ylabel(ax, ...
%     "Latitude", ...
%     "FontSize", fontSize, ...
%     "FontName", fontName);


% ================= OPTIONAL MODEL FIT =================

% % Change these values when necessary
% pearson_fit = 0.724;
% mse_fit     = 0.0343;
% 
% txt = sprintf( ...
%     "Pearson = %.3f\nMSE = %.4f", ...
%     pearson_fit, ...
%     mse_fit);
% 
% text(ax, ...
%     0.03, ...
%     0.03, ...
%     txt, ...
%     "Units", "normalized", ...
%     "HorizontalAlignment", "left", ...
%     "VerticalAlignment", "bottom", ...
%     "FontSize", 20, ...
%     "FontName", fontName);


%% ================= AXIS FORMAT =================

decorate( ...
    ax, ...
    xlimMap, ...
    ylimMap, ...
    fontSize, ...
    fontName);

hold(ax, "off");


%% ================= EXPORT =================

if ~isfolder(outdir)
    mkdir(outdir);
end


exportgraphics(fig, ...
    fullfile(outdir, "four_class_mortality_map.png"), ...
    "Resolution", 300);
savefig(fig, fullfile(outdir, "four_class_mortality_map.fig"));
exportgraphics(fig, ...
    fullfile(outdir, "four_class_mortality_map.eps"), ...
    "ContentType", "vector");

%% ================= LOCAL FUNCTION =================

function decorate(ax, xlimMap, ylimMap, fontSize, fontName)

    axis(ax, "equal");

    xlim(ax, xlimMap);
    ylim(ax, ylimMap);

    % box(ax, "on");
    % grid(ax, "on");

    ax.LineWidth = 1.0;
    ax.FontSize = fontSize;
    ax.FontName = fontName;
    ax.GridAlpha = 0.18;
    ax.TickDir = "out";

end
