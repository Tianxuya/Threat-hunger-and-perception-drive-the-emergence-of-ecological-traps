
clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================

baseDir = paths.namibia_results;

outdir = paths.namibia_figures;

if ~isfolder(outdir)
    mkdir(outdir);
end

%% ================= SCENARIOS =================

folders = {
    'uniform_low_food'
    'notoxin'
    'full'
};

labels = {
    'Uniform food (0.025737) + toxin'
    'Normal food + no toxin'
    'Normal food + toxin'
};

tags = {
    'uniform_food_toxin'
    'normal_food_no_toxin'
    'normal_food_toxin'
};

nScen = numel(folders);

%% ================= SETTINGS =================

nx = 300;
ny = 300;
distThresh = 0.12;
minAlive = 1e-15;

nClasses = 6;

fontSize = 24;
fontName = 'Times New Roman';

xlimMap = [11.99489 21.0];
ylimMap = [-29.0 -16.9059];

%% ================= BLUE COLORMAP =================

blueAnchors = [
    0.812 0.894 0.961
    0.612 0.796 0.910
    0.376 0.655 0.824
    0.157 0.490 0.698
    0.031 0.294 0.514
];

cmap = interp1( ...
    linspace(0,1,size(blueAnchors,1)), ...
    blueAnchors, ...
    linspace(0,1,nClasses), ...
    'linear');

%% ================= LOAD ALL TABLES =================

T = cell(nScen,1);
allLon = [];
allLat = [];

for s = 1:nScen

    file = fullfile(baseDir, folders{s}, ...
        'pde_mortality_grid.csv');

    if ~isfile(file)
        error('Missing file: %s', file);
    end

    T{s} = readtable(file);

    needed = {'x','y','final_hunger','final_population_mass'};

    if ~all(ismember(needed, T{s}.Properties.VariableNames))
        error('Missing required columns in %s', file);
    end

    allLon = [allLon; T{s}.x(:)]; %#ok<AGROW>
    allLat = [allLat; T{s}.y(:)]; %#ok<AGROW>

end

%% ================= COMMON INTERPOLATION GRID =================

coordValid = isfinite(allLon) & isfinite(allLat);

if ~any(coordValid)
    error('No valid coordinates found.');
end

xq = linspace( ...
    min(allLon(coordValid)), ...
    max(allLon(coordValid)), nx);

yq = linspace( ...
    min(allLat(coordValid)), ...
    max(allLat(coordValid)), ny);

[Xq, Yq] = meshgrid(xq, yq);

queryCoords = [Xq(:), Yq(:)];

%% ================= INTERPOLATE ALL SCENARIOS =================

hungerMaps = cell(nScen,1);

for s = 1:nScen

    lon = T{s}.x(:);
    lat = T{s}.y(:);
    hunger = T{s}.final_hunger(:);
    pop = T{s}.final_population_mass(:);

    valid = isfinite(lon) & isfinite(lat) & ...
            isfinite(hunger) & isfinite(pop) & ...
            (pop > minAlive);

    if nnz(valid) < 4
        error('Insufficient surviving-population cells for %s', ...
            labels{s});
    end

    lonV = lon(valid);
    latV = lat(valid);
    hungerV = hunger(valid);

    %% ---------- INTERPOLATION ----------

    F = scatteredInterpolant( ...
        lonV, latV, hungerV, ...
        'natural', 'none');

    Z = F(Xq, Yq);

    %% ---------- PRESERVE SPATIAL GAPS ----------

    idx = knnsearch( ...
        [lonV latV], queryCoords);

    dist = hypot( ...
        Xq(:) - lonV(idx), ...
        Yq(:) - latV(idx));

    nearestDist = reshape(dist, size(Xq));

    Z(nearestDist > distThresh) = NaN;

    hungerMaps{s} = Z;

end

%% ================= SHARED QUANTILE THRESHOLDS =================

% Combine interpolated hunger values across ALL scenarios

pooledValues = [];

for s = 1:nScen

    Z = hungerMaps{s};

    pooledValues = [pooledValues; ...
        Z(isfinite(Z))]; %#ok<AGROW>

end

if isempty(pooledValues)
    error('No finite interpolated hunger values.');
end

% Shared thresholds: Q25, Q50, Q75
quantileLevels = (1:nClasses-1)/nClasses;

thresholds = quantile(pooledValues, quantileLevels);

fprintf('\n=====================================\n');
fprintf('SHARED FINAL HUNGER THRESHOLDS\n');
fprintf('=====================================\n');

for k = 1:numel(thresholds)
    fprintf('Q%.0f = %.6f\n', ...
        quantileLevels(k)*100, thresholds(k));
end

if numel(unique(thresholds)) < numel(thresholds)
    warning('Repeated thresholds: fewer than 4 classes may be occupied.');
end

%% ================= PLOT THREE SCENARIOS =================

for s = 1:nScen

    Z = hungerMaps{s};

    good = isfinite(Z);

    if ~any(good(:))
        error('No finite values for %s', labels{s});
    end

    %% ---------- ASSIGN SHARED CLASSES ----------

    Zclass = nan(size(Z));

    Zclass(good) = 1;

    for k = 1:numel(thresholds)
        Zclass(good & Z > thresholds(k)) = k + 1;
    end

    %% ---------- FIGURE ----------

    fig = figure( ...
        'Color','w', ...
        'Position',[100 100 900 900], ...
        'Visible','on');

    ax = axes(fig);

    im = imagesc(ax, xq, yq, Zclass);

    set(ax, 'YDir', 'normal', 'Color', 'w');

    set(im, 'AlphaData', double(isfinite(Zclass)));

    colormap(ax, cmap);

    clim(ax, [0.5 nClasses+0.5]);

    %% ---------- SHARED COLORBAR ----------

    cb = colorbar(ax);

    cb.Ticks = 1:4;

    cb.TickLabels = { ...
        'Slight hunger', ...
        'Moderate hunger', ...
        'High hunger', ...
        'Severe hunger'};

    cb.FontSize = fontSize;
    cb.FontName = fontName;
    cb.TickLength = 0;

    cb.Label.String = 'Relative hunger level';
    cb.Label.FontName = fontName;

    %% ---------- AXES ----------

    xlabel(ax, 'Longitude', ...
        'FontSize', fontSize, ...
        'FontName', fontName);

    ylabel(ax, 'Latitude', ...
        'FontSize', fontSize, ...
        'FontName', fontName);

    title(ax, ...
        sprintf('Final spatial hunger: %s', labels{s}), ...
        'FontSize', fontSize, ...
        'FontName', fontName, ...
        'FontWeight', 'normal');

    axis(ax, 'equal');

    xlim(ax, xlimMap);
    ylim(ax, ylimMap);

    ax.LineWidth = 1.0;
    ax.FontSize = fontSize;
    ax.FontName = fontName;
    ax.TickDir = 'out';

    box(ax, 'on');

    %% ---------- EXPORT ----------

    filename = sprintf( ...
        'namibia_final_hunger_%s_%dclass_shared', ...
        tags{s}, nClasses);

    exportgraphics(fig, ...
        fullfile(outdir, [filename '.png']), ...
        'Resolution', 300);

    exportgraphics(fig, ...
        fullfile(outdir, [filename '.eps']), ...
        'ContentType', 'vector');

    savefig(fig, ...
        fullfile(outdir, [filename '.fig']));

    fprintf('Saved: %s\n', filename);

end

fprintf('\nAll three shared-class hunger maps saved in:\n%s\n', ...
    outdir);
