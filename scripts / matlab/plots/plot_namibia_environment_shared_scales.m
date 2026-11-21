clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================
baseDir = paths.namibia_results;
outdir = paths.namibia_figures;
if ~isfolder(outdir), mkdir(outdir); end

%% ================= SWITCH =================
useInterpolation = true;    % true = interpolated map
% useInterpolation = false; % false = no interpolation, original cells only

markerSize = 18;  % used only when useInterpolation = false

%% ================= THREE SCENARIOS =================
% Baseline overlap, uniform-food overlap, and food-only map.
% Food classification uses one baseline-defined green scale in ALL maps.
% Toxin uses one baseline-defined purple scale in the two overlap maps.
folders = {
    'full'
    'uniform_low_food'
    'notoxin'
};
titles = {'Food + toxin', 'Uniform food (0.025737) + toxin', 'Food only (toxin inactive)'};
tags = {'food_toxin', 'uniform_food_toxin', 'food_only'};
nScen = numel(folders);

%% ================= SETTINGS =================
fontName = 'Times New Roman';
fontSize = 18;
nx = 500; 
ny = 500;
margin = 0.02;
alphaFood = 0.50;
alphaToxin = 0.50;

cmapFood = [
    0.88 0.97 0.88
    0.66 0.87 0.66
    0.35 0.69 0.35
    0.00 0.39 0.00
];

cmapToxin = [
    0.96 0.90 0.98
    0.82 0.65 0.88
    0.62 0.35 0.72
    0.35 0.08 0.45
];

%% ================= LOAD DATA =================
T = cell(nScen,1);

for s = 1:nScen
    filepath = fullfile(baseDir, folders{s}, 'pde_mortality_grid.csv');
    if ~isfile(filepath)
        error('Missing file: %s', filepath);
    end

    T{s} = readtable(filepath);

    needed = {'x','y','food_model','toxin_model'};
    if ~all(ismember(needed, T{s}.Properties.VariableNames))
        error('Missing required columns in %s', filepath);
    end
end

% Reference thresholds use the original baseline landscape.
baselineFood = T{1}.food_model;
baselineToxin = T{1}.toxin_model;

baselineFood = baselineFood(isfinite(baselineFood));
baselineToxin = baselineToxin(isfinite(baselineToxin));

foodQ = quantile(baselineFood, [0.25 0.50 0.75]);
toxinQ = quantile(baselineToxin, [0.25 0.50 0.75]);

fprintf('Shared food quartile boundaries: %.4g %.4g %.4g\n', foodQ);
fprintf('Shared toxin quartile boundaries: %.4g %.4g %.4g\n', toxinQ);

%% ================= COMMON GEOGRAPHIC GRID =================
allLon = [];
allLat = [];

for s = 1:nScen
    allLon = [allLon; T{s}.x(:)]; %#ok<AGROW>
    allLat = [allLat; T{s}.y(:)]; %#ok<AGROW>
end

v = isfinite(allLon) & isfinite(allLat);
if ~any(v)
    error('No finite coordinates in scenario files.');
end

xq = linspace(min(allLon(v)), max(allLon(v)), nx);
yq = linspace(min(allLat(v)), max(allLat(v)), ny);
[Xq, Yq] = meshgrid(xq, yq);
xyQuery = [Xq(:) Yq(:)];

xlims = [xq(1) xq(end)] + [-1 1] * margin * (xq(end) - xq(1));
ylims = [yq(1) yq(end)] + [-1 1] * margin * (yq(end) - yq(1));

%% ================= PREPARE MAP DATA =================
RGB = cell(nScen,1);
validMap = cell(nScen,1);

cellLon = cell(nScen,1);
cellLat = cell(nScen,1);
cellRGB = cell(nScen,1);

for s = 1:nScen

    lon = T{s}.x(:);
    lat = T{s}.y(:);
    food = T{s}.food_model(:);
    toxin = T{s}.toxin_model(:);

    valid = isfinite(lon) & isfinite(lat) & isfinite(food);
    if s < 3
        valid = valid & isfinite(toxin);
    end

    lon = lon(valid);
    lat = lat(valid);
    food = food(valid);
    toxin = toxin(valid);

    if numel(food) < 4
        error('Insufficient cells in %s', titles{s});
    end

    %% ---------- classify original cell values ----------
    Sclass_cell = ones(size(food));
    for k = 1:3
        Sclass_cell(food > foodQ(k)) = k + 1;
    end

    if s < 3
        Tclass_cell = ones(size(toxin));
        for k = 1:3
            Tclass_cell(toxin > toxinQ(k)) = k + 1;
        end
    else
        Tclass_cell = nan(size(food));
    end

    %% ---------- build RGB at original cells ----------
    thisRGB = zeros(numel(food), 3);

    for channel = 1:3
        foodChannel = cmapFood(Sclass_cell, channel);

        if s < 3
            toxinChannel = cmapToxin(Tclass_cell, channel);
            thisRGB(:,channel) = alphaFood * foodChannel + alphaToxin * toxinChannel;
        else
            thisRGB(:,channel) = foodChannel;
        end
    end

    cellLon{s} = lon;
    cellLat{s} = lat;
    cellRGB{s} = thisRGB;

    %% ---------- interpolated version ----------
    if useInterpolation

        Fs = scatteredInterpolant(lon, lat, food, 'natural', 'none');
        Zfood = Fs(Xq, Yq);

        if s < 3
            Ft = scatteredInterpolant(lon, lat, toxin, 'natural', 'none');
            Ztoxin = Ft(Xq, Yq);
        else
            Ztoxin = nan(size(Xq));
        end

        % Same original no-data masking rule
        [~,neighborDist] = knnsearch([lon lat], [lon lat], 'K', 2);
        cellSpacing = median(neighborDist(:,2), 'omitnan');

        idx = knnsearch([lon lat], xyQuery);
        nearestDist = reshape(hypot(Xq(:)-lon(idx), Yq(:)-lat(idx)), size(Xq));

        mask = nearestDist <= 1.8 * cellSpacing & isfinite(Zfood);
        if s < 3
            mask = mask & isfinite(Ztoxin);
        end

        Sclass = ones(size(Zfood));
        for k = 1:3
            Sclass(Zfood > foodQ(k)) = k + 1;
        end

        Tclass = ones(size(Ztoxin));
        if s < 3
            for k = 1:3
                Tclass(Ztoxin > toxinQ(k)) = k + 1;
            end
        end

        imRGB = ones([size(Xq) 3]);

        for channel = 1:3
            foodChannel = reshape(cmapFood(Sclass, channel), size(Xq));

            if s < 3
                toxinChannel = reshape(cmapToxin(Tclass, channel), size(Xq));
                imRGB(:,:,channel) = alphaFood * foodChannel + alphaToxin * toxinChannel;
            else
                imRGB(:,:,channel) = foodChannel;
            end
        end

        RGB{s} = imRGB;
        validMap{s} = mask;
    end
end

%% ================= OVERLAP KEY =================
keyRGB = zeros(4,4,3);

for i = 1:4
    for j = 1:4
        keyRGB(i,j,:) = ...
            alphaFood  * reshape(cmapFood(j,:),1,1,3) + ...
            alphaToxin * reshape(cmapToxin(i,:),1,1,3);
    end
end

%% ================= THREE SEPARATE FIGURES =================
for s = 1:nScen

    fig = figure('Color','w','Position',[100 100 1000 850],'Visible','on');
    ax = axes(fig);
    hold(ax, 'on');

    if useInterpolation
        him = image(ax, [xq(1) xq(end)], [yq(1) yq(end)], RGB{s});
        set(ax, 'YDir', 'normal', 'Color', 'w');
        set(him, 'AlphaData', double(validMap{s}));
    else
        scatter(ax, cellLon{s}, cellLat{s}, markerSize, cellRGB{s}, ...
            's', 'filled', 'MarkerEdgeColor', 'none');
        set(ax, 'YDir', 'normal', 'Color', 'w');
    end

    axis(ax, 'equal');
    xlim(ax, xlims);
    ylim(ax, ylims);

    xlabel(ax, 'Longitude', 'FontName', fontName, 'FontSize', fontSize);
    ylabel(ax, 'Latitude',  'FontName', fontName, 'FontSize', fontSize);
    title(ax, titles{s}, 'FontName', fontName, 'FontSize', fontSize, 'FontWeight', 'normal');

    ax.LineWidth = 1;
    ax.FontSize = fontSize;
    ax.FontName = fontName;
    ax.GridAlpha = 0.15;
    ax.TickDir = 'out';

    box(ax, 'on');
    grid(ax, 'on');

    % Food colorbar (conceptual shared scale)
    colormap(ax, cmapFood);
    clim(ax, [0.5 4.5]);

    cb = colorbar(ax);
    cb.Ticks = 1:4;
    cb.TickLabels = {'Very low','Low','Medium','High'};
    cb.FontName = fontName;
    cb.FontSize = fontSize;
    cb.TickLength = 0;
    cb.Label.String = 'Food level (shared baseline classes)';
    cb.Label.FontName = fontName;

    if s < 3
        % Shared 4x4 food-toxin key
        axKey = axes(fig, 'Position', [.19 .20 .19 .19]);
        image(axKey, 1:4, 1:4, keyRGB);
        set(axKey, 'YDir', 'normal', ...
            'XTick', [1 4], 'YTick', [1 4], ...
            'XTickLabel', {'L','H'}, ...
            'YTickLabel', {'L','H'}, ...
            'FontSize', 15, ...
            'FontName', fontName, ...
            'LineWidth', 1);

        xlabel(axKey, 'Food',  'FontName', fontName, 'FontSize', 15);
        ylabel(axKey, 'Toxin', 'FontName', fontName, 'FontSize', 15);
        box(axKey, 'on');
    end

    if useInterpolation
        basename = sprintf('namibia_environment_%s_interpolated', tags{s});
    else
        basename = sprintf('namibia_environment_%s_no_interpolation', tags{s});
    end

    exportgraphics(fig, fullfile(outdir, [basename '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(outdir, [basename '.eps']), 'ContentType', 'vector');
    savefig(fig, fullfile(outdir, [basename '.fig']));

    fprintf('Saved %s (.png, .eps, .fig)\n', basename);
end

fprintf('All three maps saved to: %s\n', outdir);
