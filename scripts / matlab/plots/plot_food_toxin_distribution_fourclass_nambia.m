clear; clc; close all;
paths = project_paths();

%% ================================================================
% FOUR-CLASS MAPS FOR
%
%   s(x)   = food_model
%   tau(x) = toxin_model
%
% PLUS:
%   Food + toxin overlap map
%
% Food  = green
% Toxin = purple
%
% All figures are generated invisibly.
%% ================================================================

%% ================= PATHS =================

indir = fullfile(paths.namibia_results, 'full');
outdir = paths.namibia_figures;
if ~isfolder(outdir), mkdir(outdir); end

gridFile = fullfile(indir, "pde_mortality_grid.csv");

if ~isfile(gridFile)
    error("Cannot find:\n%s", gridFile);
end

%% ================= READ GRID =================

PDEgrid = readtable(gridFile);

fprintf("\nAvailable variables:\n");
disp(PDEgrid.Properties.VariableNames');

%% ================= COORDINATES =================

lon = PDEgrid.x;
lat = PDEgrid.y;

%% ================= FIGURE SETTINGS =================

fontSize = 50;
fontName = "Arial";

nx = 500;
ny = 500;

margin = 0.02;

%% ================================================================
%% FOOD s(x)
%% ================================================================

if ~ismember("food_model", PDEgrid.Properties.VariableNames)
    error("Column food_model does not exist in pde_mortality_grid.csv");
end

s = PDEgrid.food_model;

%% ================= CLEAN FOOD DATA =================

valid_s = ...
    isfinite(lon) & ...
    isfinite(lat) & ...
    isfinite(s);

lon_s = lon(valid_s);
lat_s = lat(valid_s);
s     = s(valid_s);

%% ================= FOUR FOOD CLASSES =================

q1_s = quantile(s, 0.25);
q2_s = quantile(s, 0.50);
q3_s = quantile(s, 0.75);

fprintf("\n========================================\n");
fprintf("FOOD s(x) FOUR-CLASS THRESHOLDS\n");
fprintf("========================================\n");

fprintf("Q1 = %.6f\n", q1_s);
fprintf("Q2 = %.6f\n", q2_s);
fprintf("Q3 = %.6f\n", q3_s);

%% ================= FOOD INTERPOLATION GRID =================

xq_s = linspace(min(lon_s), max(lon_s), nx);
yq_s = linspace(min(lat_s), max(lat_s), ny);

[Xq_s, Yq_s] = meshgrid(xq_s, yq_s);

%% ================= FOOD INTERPOLATION =================

F_s = scatteredInterpolant( ...
    lon_s, ...
    lat_s, ...
    s, ...
    "natural", ...
    "none");

Zs = F_s(Xq_s, Yq_s);

%% ================= PRESERVE NO-DATA REGIONS =================

[~, modelNeighborDist_s] = knnsearch( ...
    [lon_s lat_s], ...
    [lon_s lat_s], ...
    "K", 2);

cellSpacing_s = median( ...
    modelNeighborDist_s(:,2), ...
    "omitnan");

idx_s = knnsearch( ...
    [lon_s lat_s], ...
    [Xq_s(:) Yq_s(:)]);

nearestDist_s = sqrt( ...
    (Xq_s(:) - lon_s(idx_s)).^2 + ...
    (Yq_s(:) - lat_s(idx_s)).^2);

nearestDist_s = reshape( ...
    nearestDist_s, ...
    size(Xq_s));

distThresh_s = 1.8 * cellSpacing_s;

Zs(nearestDist_s > distThresh_s) = NaN;

%% ================= CLASSIFY FOOD =================

validZs = isfinite(Zs);

Sclass = nan(size(Zs));

% Class 1
Sclass( ...
    validZs & ...
    Zs <= q1_s) = 1;

% Class 2
Sclass( ...
    validZs & ...
    Zs > q1_s & ...
    Zs <= q2_s) = 2;

% Class 3
Sclass( ...
    validZs & ...
    Zs > q2_s & ...
    Zs <= q3_s) = 3;

% Class 4
Sclass( ...
    validZs & ...
    Zs > q3_s) = 4;

%% ================= FOOD COLORMAP =================

% Green: light -> dark
cmap_s = [
    0.88 0.97 0.88   % Very low
    0.66 0.87 0.66   % Low
    0.35 0.69 0.35   % Medium
    0.00 0.39 0.00   % High
];

%% ================= PLOT FOOD =================

figS = figure( ...
    "Color", "w", ...
    "Position", [100 100 1000 850], ...
    "Visible", "off");

axS = axes(figS);

contourf( ...
    axS, ...
    Xq_s, ...
    Yq_s, ...
    Sclass, ...
    [1 2 3 4 5], ...
    "LineColor", "none");

axS.Color = "w";

colormap(axS, cmap_s);

clim(axS, [1 5]);

%% ================= FOOD COLORBAR =================

cbS = colorbar(axS);

cbS.Ticks = [1.5 2.5 3.5 4.5];

cbS.TickLabels = { ...
    "Very low", ...
    "Low", ...
    "Medium", ...
    "High"};

cbS.FontSize = fontSize;
cbS.FontName = fontName;
cbS.TickLength = 0;

%% ================= FOOD AXES =================

xmin_s = min(lon_s);
xmax_s = max(lon_s);

ymin_s = min(lat_s);
ymax_s = max(lat_s);

dx_s = xmax_s - xmin_s;
dy_s = ymax_s - ymin_s;

xlim(axS, [ ...
    xmin_s - margin*dx_s, ...
    xmax_s + margin*dx_s]);

ylim(axS, [ ...
    ymin_s - margin*dy_s, ...
    ymax_s + margin*dy_s]);

axis(axS, "equal");

xlabel(axS, ...
    "Longitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

ylabel(axS, ...
    "Latitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

box(axS, "on");
grid(axS, "on");

axS.LineWidth = 1.0;
axS.FontSize = fontSize;
axS.FontName = fontName;
axS.GridAlpha = 0.15;
axS.TickDir = "out";

%% ================= EXPORT FOOD =================

exportgraphics(figS, fullfile(outdir, "four_class_food_s_map.png"), ...
    "Resolution", 300);
exportgraphics(figS, fullfile(outdir, "four_class_food_s_map.eps"), ...
    "ContentType", "vector");
savefig(figS, fullfile(outdir, "four_class_food_s_map.fig"));

%% ================================================================
%% TOXIN tau(x)
%% ================================================================

if ~ismember("toxin_model", PDEgrid.Properties.VariableNames)
    error("Column toxin_model does not exist in pde_mortality_grid.csv");
end

tau = PDEgrid.toxin_model;

%% ================= CLEAN TOXIN DATA =================

valid_t = ...
    isfinite(lon) & ...
    isfinite(lat) & ...
    isfinite(tau);

lon_t = lon(valid_t);
lat_t = lat(valid_t);
tau   = tau(valid_t);

%% ================= FOUR TOXIN CLASSES =================

q1_t = quantile(tau, 0.25);
q2_t = quantile(tau, 0.50);
q3_t = quantile(tau, 0.75);

fprintf("\n========================================\n");
fprintf("TOXIN tau(x) FOUR-CLASS THRESHOLDS\n");
fprintf("========================================\n");

fprintf("Q1 = %.6f\n", q1_t);
fprintf("Q2 = %.6f\n", q2_t);
fprintf("Q3 = %.6f\n", q3_t);

%% ================= TOXIN INTERPOLATION GRID =================

xq_t = linspace(min(lon_t), max(lon_t), nx);
yq_t = linspace(min(lat_t), max(lat_t), ny);

[Xq_t, Yq_t] = meshgrid(xq_t, yq_t);

%% ================= TOXIN INTERPOLATION =================

F_t = scatteredInterpolant( ...
    lon_t, ...
    lat_t, ...
    tau, ...
    "natural", ...
    "none");

Zt = F_t(Xq_t, Yq_t);

%% ================= PRESERVE NO-DATA REGIONS =================

[~, modelNeighborDist_t] = knnsearch( ...
    [lon_t lat_t], ...
    [lon_t lat_t], ...
    "K", 2);

cellSpacing_t = median( ...
    modelNeighborDist_t(:,2), ...
    "omitnan");

idx_t = knnsearch( ...
    [lon_t lat_t], ...
    [Xq_t(:) Yq_t(:)]);

nearestDist_t = sqrt( ...
    (Xq_t(:) - lon_t(idx_t)).^2 + ...
    (Yq_t(:) - lat_t(idx_t)).^2);

nearestDist_t = reshape( ...
    nearestDist_t, ...
    size(Xq_t));

distThresh_t = 1.8 * cellSpacing_t;

Zt(nearestDist_t > distThresh_t) = NaN;

%% ================= CLASSIFY TOXIN =================

validZt = isfinite(Zt);

Tclass = nan(size(Zt));

% Class 1
Tclass( ...
    validZt & ...
    Zt <= q1_t) = 1;

% Class 2
Tclass( ...
    validZt & ...
    Zt > q1_t & ...
    Zt <= q2_t) = 2;

% Class 3
Tclass( ...
    validZt & ...
    Zt > q2_t & ...
    Zt <= q3_t) = 3;

% Class 4
Tclass( ...
    validZt & ...
    Zt > q3_t) = 4;

%% ================= TOXIN COLORMAP =================

% Purple: light -> dark
cmap_t = [
    0.96 0.90 0.98   % Very low
    0.82 0.65 0.88   % Low
    0.62 0.35 0.72   % Medium
    0.35 0.08 0.45   % High
];

%% ================= PLOT TOXIN =================

figT = figure( ...
    "Color", "w", ...
    "Position", [150 100 1000 850], ...
    "Visible", "off");

axT = axes(figT);

contourf( ...
    axT, ...
    Xq_t, ...
    Yq_t, ...
    Tclass, ...
    [1 2 3 4 5], ...
    "LineColor", "none");

axT.Color = "w";

colormap(axT, cmap_t);

clim(axT, [1 5]);

%% ================= TOXIN COLORBAR =================

cbT = colorbar(axT);

cbT.Ticks = [1.5 2.5 3.5 4.5];

cbT.TickLabels = { ...
    "Very low", ...
    "Low", ...
    "Medium", ...
    "High"};

cbT.FontSize = fontSize;
cbT.FontName = fontName;
cbT.TickLength = 0;

%% ================= TOXIN AXES =================

xmin_t = min(lon_t);
xmax_t = max(lon_t);

ymin_t = min(lat_t);
ymax_t = max(lat_t);

dx_t = xmax_t - xmin_t;
dy_t = ymax_t - ymin_t;

xlim(axT, [ ...
    xmin_t - margin*dx_t, ...
    xmax_t + margin*dx_t]);

ylim(axT, [ ...
    ymin_t - margin*dy_t, ...
    ymax_t + margin*dy_t]);

axis(axT, "equal");

xlabel(axT, ...
    "Longitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

ylabel(axT, ...
    "Latitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

box(axT, "on");
grid(axT, "on");

axT.LineWidth = 1.0;
axT.FontSize = fontSize;
axT.FontName = fontName;
axT.GridAlpha = 0.15;
axT.TickDir = "out";

%% ================= EXPORT TOXIN =================

exportgraphics( ...
    figT, ...
    fullfile(outdir, "four_class_toxin_tau_map.png"), ...
    "Resolution", 300);

exportgraphics( ...
    figT, ...
    fullfile(outdir, "four_class_toxin_tau_map.eps"), ...
    "ContentType", "vector");

savefig( ...
    figT, ...
    fullfile(outdir, "four_class_toxin_tau_map.fig"));

%% ================================================================
%% FOOD + TOXIN OVERLAP MAP
%% ================================================================

valid_both = ...
    isfinite(lon) & ...
    isfinite(lat) & ...
    isfinite(PDEgrid.food_model) & ...
    isfinite(PDEgrid.toxin_model);

lon_b = lon(valid_both);
lat_b = lat(valid_both);

s_b   = PDEgrid.food_model(valid_both);
tau_b = PDEgrid.toxin_model(valid_both);

%% ================= COMMON INTERPOLATION GRID =================

xq = linspace(min(lon_b), max(lon_b), nx);
yq = linspace(min(lat_b), max(lat_b), ny);

[Xq, Yq] = meshgrid(xq, yq);

%% ================= INTERPOLATE FOOD =================

Fs = scatteredInterpolant( ...
    lon_b, ...
    lat_b, ...
    s_b, ...
    "natural", ...
    "none");

Zs_overlap = Fs(Xq, Yq);

%% ================= INTERPOLATE TOXIN =================

Ft = scatteredInterpolant( ...
    lon_b, ...
    lat_b, ...
    tau_b, ...
    "natural", ...
    "none");

Zt_overlap = Ft(Xq, Yq);

%% ================= PRESERVE DOMAIN / HOLES =================

[~, modelNeighborDist_b] = knnsearch( ...
    [lon_b lat_b], ...
    [lon_b lat_b], ...
    "K", 2);

cellSpacing_b = median( ...
    modelNeighborDist_b(:,2), ...
    "omitnan");

idx_b = knnsearch( ...
    [lon_b lat_b], ...
    [Xq(:) Yq(:)]);

nearestDist_b = sqrt( ...
    (Xq(:) - lon_b(idx_b)).^2 + ...
    (Yq(:) - lat_b(idx_b)).^2);

nearestDist_b = reshape( ...
    nearestDist_b, ...
    size(Xq));

distThresh_b = 1.8 * cellSpacing_b;

Zs_overlap(nearestDist_b > distThresh_b) = NaN;
Zt_overlap(nearestDist_b > distThresh_b) = NaN;

%% ================= CLASSIFY OVERLAP FOOD =================

q1_so = quantile(s_b, 0.25);
q2_so = quantile(s_b, 0.50);
q3_so = quantile(s_b, 0.75);

validS = isfinite(Zs_overlap);

Sclass_overlap = nan(size(Zs_overlap));

Sclass_overlap( ...
    validS & ...
    Zs_overlap <= q1_so) = 1;

Sclass_overlap( ...
    validS & ...
    Zs_overlap > q1_so & ...
    Zs_overlap <= q2_so) = 2;

Sclass_overlap( ...
    validS & ...
    Zs_overlap > q2_so & ...
    Zs_overlap <= q3_so) = 3;

Sclass_overlap( ...
    validS & ...
    Zs_overlap > q3_so) = 4;

%% ================= CLASSIFY OVERLAP TOXIN =================

q1_to = quantile(tau_b, 0.25);
q2_to = quantile(tau_b, 0.50);
q3_to = quantile(tau_b, 0.75);

validT = isfinite(Zt_overlap);

Tclass_overlap = nan(size(Zt_overlap));

Tclass_overlap( ...
    validT & ...
    Zt_overlap <= q1_to) = 1;

Tclass_overlap( ...
    validT & ...
    Zt_overlap > q1_to & ...
    Zt_overlap <= q2_to) = 2;

Tclass_overlap( ...
    validT & ...
    Zt_overlap > q2_to & ...
    Zt_overlap <= q3_to) = 3;

Tclass_overlap( ...
    validT & ...
    Zt_overlap > q3_to) = 4;

%% ================================================================
%% BUILD FOOD RGB IMAGE
%% ================================================================

RGB_food = ones([size(Sclass_overlap), 3]);

for k = 1:4

    mask = Sclass_overlap == k;

    for c = 1:3

        temp = RGB_food(:,:,c);

        temp(mask) = cmap_s(k,c);

        RGB_food(:,:,c) = temp;

    end

end

%% ================================================================
%% BUILD TOXIN RGB IMAGE
%% ================================================================

RGB_toxin = ones([size(Tclass_overlap), 3]);

for k = 1:4

    mask = Tclass_overlap == k;

    for c = 1:3

        temp = RGB_toxin(:,:,c);

        temp(mask) = cmap_t(k,c);

        RGB_toxin(:,:,c) = temp;

    end

end

%% ================= BLEND FOOD + TOXIN =================

alphaFood  = 0.50;
alphaToxin = 0.50;

RGB_overlap = ...
    alphaFood  * RGB_food + ...
    alphaToxin * RGB_toxin;

RGB_overlap = max(0, min(1, RGB_overlap));

%% ================= VALID OVERLAP REGION =================

validOverlap = ...
    isfinite(Sclass_overlap) & ...
    isfinite(Tclass_overlap);

%% ================= PLOT OVERLAP =================

figO = figure( ...
    "Color", "w", ...
    "Position", [100 100 1000 850], ...
    "Visible", "off");

axO = axes(figO);

h = image( ...
    axO, ...
    [xq(1) xq(end)], ...
    [yq(1) yq(end)], ...
    RGB_overlap);

set(axO, "YDir", "normal");

% Hide no-data area
set(h, ...
    "AlphaData", ...
    double(validOverlap));

axO.Color = "w";

hold(axO, "on");

%% ================= CORRECT OVERLAP AXES =================

xmin_o = min(lon_b);
xmax_o = max(lon_b);

ymin_o = min(lat_b);
ymax_o = max(lat_b);

dx_o = xmax_o - xmin_o;
dy_o = ymax_o - ymin_o;

xlim(axO, [ ...
    xmin_o - margin*dx_o, ...
    xmax_o + margin*dx_o]);

ylim(axO, [ ...
    ymin_o - margin*dy_o, ...
    ymax_o + margin*dy_o]);

axis(axO, "equal");

xlabel(axO, ...
    "Longitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

ylabel(axO, ...
    "Latitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

box(axO, "on");
grid(axO, "on");

axO.LineWidth = 1.0;
axO.FontSize = fontSize;
axO.FontName = fontName;
axO.GridAlpha = 0.15;
axO.TickDir = "out";

%% ================================================================
%% 2D FOOD--TOXIN OVERLAP KEY
%% ================================================================

nFood  = 4;
nToxin = 4;

keyRGB = zeros(nToxin, nFood, 3);

for i = 1:nToxin

    for j = 1:nFood

        thisFood = reshape( ...
            cmap_s(j,:), ...
            1, 1, 3);

        thisToxin = reshape( ...
            cmap_t(i,:), ...
            1, 1, 3);

        keyRGB(i,j,:) = ...
            alphaFood  * thisFood + ...
            alphaToxin * thisToxin;

    end

end

%% ================= CREATE OVERLAP KEY =================

axKey = axes( ...
    figO, ...
    "Position", [0.20 0.20 0.22 0.22]);

image( ...
    axKey, ...
    1:nFood, ...
    1:nToxin, ...
    keyRGB);

set(axKey, "YDir", "normal");

%% ================= ONLY LOW / HIGH LABELS =================

axKey.XTick = [1 4];
axKey.YTick = [1 4];

axKey.XTickLabel = {"L","H"};
axKey.YTickLabel = {"L","H"};

axKey.FontSize = 50;
axKey.FontName = fontName;
axKey.LineWidth = 1.0;

box(axKey, "on");

xlabel( ...
    axKey, ...
    "Food", ...
    "FontSize", 50, ...
    "FontName", fontName, ...
    "FontWeight", "normal");

ylabel( ...
    axKey, ...
    "Toxin", ...
    "FontSize", 50, ...
    "FontName", fontName, ...
    "FontWeight", "normal");

hold(axO, "off");

%% ================= EXPORT OVERLAP =================

exportgraphics( ...
    figO, ...
    fullfile(outdir, "food_toxin_overlap_map.png"), ...
    "Resolution", 300);

exportgraphics( ...
    figO, ...
    fullfile(outdir, "food_toxin_overlap_map.eps"), ...
    "ContentType", "vector");

savefig( ...
    figO, ...
    fullfile(outdir, "food_toxin_overlap_map.fig"));

%% ================= FINISHED =================

fprintf("\n========================================\n");
fprintf("ALL FIGURES SAVED\n");
fprintf("========================================\n\n");

fprintf("Food map:\n");
fprintf("  four_class_food_s_map.png\n");
fprintf("  four_class_food_s_map.eps\n");
fprintf("  four_class_food_s_map.fig\n\n");

fprintf("Toxin map:\n");
fprintf("  four_class_toxin_tau_map.png\n");
fprintf("  four_class_toxin_tau_map.eps\n");
fprintf("  four_class_toxin_tau_map.fig\n\n");

fprintf("Food + toxin overlap:\n");
fprintf("  food_toxin_overlap_map.png\n");
fprintf("  food_toxin_overlap_map.eps\n");
fprintf("  food_toxin_overlap_map.fig\n\n");

fprintf("Output directory:\n%s\n\n", outdir);

%% ================= CLOSE INVISIBLE FIGURES =================

close(figS);
close(figT);
close(figO);
