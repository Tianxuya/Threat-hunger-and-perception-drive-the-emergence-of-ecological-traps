clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================

indir = fullfile(paths.namibia_results, 'full');
outdir = paths.namibia_figures;

%% ================= READ DATA =================

PDEgrid = readtable(fullfile(indir, "pde_mortality_grid.csv"));

lon = PDEgrid.x;
lat = PDEgrid.y;
initial = PDEgrid.initial_population_mass;

%% ================= CLEAN DATA =================

valid = isfinite(lon) & isfinite(lat) & isfinite(initial);

lon = lon(valid);
lat = lat(valid);
initial = initial(valid);

%% ================= OPTIONAL CLIP EXTREMES =================

vmax = prctile(initial, 99);
initial_clip = min(initial, vmax);

%% ================= INTERPOLATION GRID =================

nx = 300;
ny = 300;

xq = linspace(min(lon), max(lon), nx);
yq = linspace(min(lat), max(lat), ny);

[Xq, Yq] = meshgrid(xq, yq);

%% ================= INTERPOLATE =================

F = scatteredInterpolant( ...
    lon, lat, initial_clip, ...
    "natural", "none");

Zq = F(Xq, Yq);

%% ================= PRESERVE NO-DATA REGIONS =================

idx = knnsearch([lon lat], [Xq(:) Yq(:)]);

nearestDist = sqrt( ...
    (Xq(:) - lon(idx)).^2 + ...
    (Yq(:) - lat(idx)).^2);

nearestDist = reshape(nearestDist, size(Xq));

distThresh = 0.12;

Zq(nearestDist > distThresh) = NaN;

%% ================= CREATE TEN CLASSES =================

validZ = isfinite(Zq);
vals = Zq(validZ);

% Actual minimum and maximum of displayed data
realMin = min(vals);
realMax = max(vals);

% Decile thresholds
q = quantile(vals, 0.1:0.1:0.9);

fprintf("\nInitial-condition ten-class thresholds:\n");

for k = 1:9
    fprintf("Q%d = %.6g\n", k, q(k));
end

fprintf("\nDisplayed minimum = %.6g\n", realMin);
fprintf("Displayed maximum = %.6g\n\n", realMax);

%% ================= CLASSIFY =================

Zclass = nan(size(Zq));

% Class 1
Zclass(validZ & Zq <= q(1)) = 1;

% Classes 2--9
for k = 2:9

    Zclass( ...
        validZ & ...
        Zq > q(k-1) & ...
        Zq <= q(k)) = k;

end

% Class 10
Zclass(validZ & Zq > q(9)) = 10;

%% ================= CHECK CLASS COUNTS =================

fprintf("Class distribution:\n");

for k = 1:10

    n_k = sum(Zclass(:) == k);

    fprintf("Class %2d: %d pixels\n", k, n_k);

end

%% ================= FIGURE SETTINGS =================

fontSize = 18;
fontName = "Arial";

xlimMap = [11.99489 21.0];
ylimMap = [-29.0 -16.9059];

fig = figure( ...
    "Color", "w", ...
    "Position", [100 100 900 900]);

ax = axes(fig);

%% ================= PLOT TEN CLASSES =================

contourf(ax, ...
    Xq, Yq, Zclass, ...
    0.5:1:10.5, ...
    "LineColor", "none");

hold(ax, "on");

%% ================= TEN-COLOR COLORMAP =================

% Light yellow -> orange -> red
cmap = interp1( ...
    [1 5 10], ...
    [1.00 1.00 0.78; ...
     1.00 0.60 0.25; ...
     0.70 0.05 0.05], ...
    1:10);

colormap(ax, cmap);

clim(ax, [0.5 10.5]);

%% ================= COLORBAR =================

cb = colorbar(ax);

% Only show minimum and maximum
cb.Ticks = [1 10];

cb.TickLabels = { ...
    sprintf('%.3g', realMin), ...
    sprintf('%.3g', realMax) ...
};

cb.FontSize = fontSize;
cb.FontName = fontName;
cb.TickLength = 0;

%% ================= AXIS LABELS =================

xlabel(ax, ...
    "Longitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

ylabel(ax, ...
    "Latitude", ...
    "FontSize", fontSize, ...
    "FontName", fontName);

%% ================= AXIS FORMAT =================

decorate(ax, xlimMap, ylimMap, fontSize, fontName);

hold(ax, "off");

%% ================= EXPORT =================

if ~isfolder(outdir)
    mkdir(outdir);
end

exportgraphics(fig, ...
    fullfile(outdir, "initial_population.eps"), ...
    "ContentType", "vector");

exportgraphics(fig, ...
    fullfile(outdir, "initial_population.png"), ...
    "Resolution", 300);

savefig(fig, ...
    fullfile(outdir, "initial_population.fig"));

%% ================= LOCAL FUNCTION =================

function decorate(ax, xlimMap, ylimMap, fontSize, fontName)

    axis(ax, "equal");

    xlim(ax, xlimMap);
    ylim(ax, ylimMap);

    box(ax, "on");
    grid(ax, "on");

    ax.LineWidth = 1.0;
    ax.FontSize = fontSize;
    ax.FontName = fontName;
    ax.GridAlpha = 0.18;
    ax.TickDir = "out";

end