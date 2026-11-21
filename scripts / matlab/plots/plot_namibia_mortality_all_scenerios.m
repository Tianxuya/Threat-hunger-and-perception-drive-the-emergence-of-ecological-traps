clear; clc; close all;
paths = project_paths();

%% ================= PATHS =================
baseDir = paths.namibia_results;
outdir = paths.namibia_figures;
if ~isfolder(outdir), mkdir(outdir); end

%% ================= INTERPOLATION SWITCH =================
useInterpolation = true;  % true: natural-neighbor interpolation; false: original spatial cells
markerSize = 18;          % square marker area in points^2, used only when false

%% ================= THREE SCENARIOS =================
folders = {
    'uniform_low_food'
    'notoxin'
    'full'
};
tags = {'uniform_food_toxin','normal_food_no_toxin','normal_food_toxin'};
titles = {'Uniform mean food (0.025737) + toxin', ...
    'Normal food + no toxin', ...
    'Normal food + toxin'};

%% ================= SETTINGS =================
fontSize = 24;
fontName = 'Times New Roman';
xlimMap = [11.99489 21.0];
ylimMap = [-29.0 -16.9059];
nx = 300;
ny = 300;
distThresh = 0.12;
showModelBackground = true;
showBubbleLabels = true;

cmap = [
    1.00 1.00 0.75
    1.00 0.85 0.45
    1.00 0.60 0.30
    0.80 0.10 0.10
];
sizeClass = [35 55 80 110 150 200 270 350];
obsLabels = {'1-3','4-6','7-10','11-16', ...
    '17-24','25-33','34-50','51-79'};

%% ================= READ MODEL AND OBSERVATIONS =================
nScen = numel(folders);
modelMaps = cell(nScen,1);
mortTable = cell(nScen,1);
lonAll = cell(nScen,1);
latAll = cell(nScen,1);
xqAll = cell(nScen,1);
yqAll = cell(nScen,1);
rawDeath = cell(nScen,1);
allDeathValues = [];

for i = 1:nScen
    indir = fullfile(baseDir,folders{i});
    gridFile = fullfile(indir,'pde_mortality_grid.csv');
    mortFile = fullfile(indir,'mortality_event_comparison.csv');
    if ~isfile(gridFile), error('Missing file: %s',gridFile); end
    if ~isfile(mortFile), error('Missing file: %s',mortFile); end

    PDEgrid = readtable(gridFile);
    mort = readtable(mortFile);
    neededGrid = {'x','y','cumulative_total_mortality_mass'};
    neededObs = {'X','Y','N_Dead'};
    if ~all(ismember(neededGrid,PDEgrid.Properties.VariableNames))
        error('Missing model data columns in %s',gridFile);
    end
    if ~all(ismember(neededObs,mort.Properties.VariableNames))
        error('Missing observed data columns in %s',mortFile);
    end

    lon = PDEgrid.x(:);
    lat = PDEgrid.y(:);
    death = PDEgrid.cumulative_total_mortality_mass(:);
    good = isfinite(lon) & isfinite(lat) & isfinite(death);
    lon = lon(good); lat = lat(good); death = death(good);
    if numel(death) < 4, error('Not enough model points in %s',gridFile); end
    lonAll{i} = lon;
    latAll{i} = lat;
    rawDeath{i} = death;
    allDeathValues = [allDeathValues;death]; %#ok<AGROW>

    goodObs = isfinite(mort.X) & isfinite(mort.Y) & isfinite(mort.N_Dead) & ...
        mort.N_Dead >= 1 & mort.N_Dead < 80;
    mortTable{i} = mort(goodObs,:);
end

%% ================= SHARED 99% CUTOFF =================
sharedVmax = prctile(allDeathValues,99);
fprintf('Shared 99%% mortality cutoff = %.6g\n',sharedVmax);

%% ================= PREPARE MAPS =================
pooled = [];
for i = 1:nScen
    lon = lonAll{i};
    lat = latAll{i};
    death_clip = min(rawDeath{i},sharedVmax);

    if useInterpolation
        xq = linspace(min(lon),max(lon),nx);
        yq = linspace(min(lat),max(lat),ny);
        [Xq,Yq] = meshgrid(xq,yq);
        F = scatteredInterpolant(lon,lat,death_clip,'natural','none');
        Z = F(Xq,Yq);
        idx = knnsearch([lon lat],[Xq(:) Yq(:)]);
        nearestDist = reshape(hypot(Xq(:)-lon(idx),Yq(:)-lat(idx)),size(Xq));
        Z(nearestDist > distThresh) = NaN;
        xqAll{i} = xq;
        yqAll{i} = yq;
    else
        Z = death_clip;  % Original node values, without interpolation
    end

    modelMaps{i} = Z;
    pooled = [pooled; Z(isfinite(Z))]; %#ok<AGROW>
end

%% ================= SHARED FOUR MODEL CLASSES =================
if isempty(pooled), error('No finite modeled mortality values.'); end
q = quantile(pooled,[0.25 0.50 0.75]);
fprintf('Shared model mortality quartiles: %.6g, %.6g, %.6g\n',q);

%% ================= PLOT AND SAVE ALL SCENARIOS =================
for i = 1:nScen
    fig = figure('Color','w','Position',[100 100 900 900],'Visible','on');
    ax = axes(fig);
    hold(ax,'on');

    if showModelBackground
        Z = modelMaps{i};
        validZ = isfinite(Z);
        Zclass = nan(size(Z));
        Zclass(validZ) = 1;
        Zclass(validZ & Z > q(1)) = 2;
        Zclass(validZ & Z > q(2)) = 3;
        Zclass(validZ & Z > q(3)) = 4;

        if useInterpolation
            h = imagesc(ax,xqAll{i},yqAll{i},Zclass);
            set(h,'AlphaData',double(validZ));
        else
            scatter(ax,lonAll{i},latAll{i},markerSize,Zclass, ...
                's','filled','MarkerEdgeColor','none');
        end
        set(ax,'YDir','normal','Color','w');
        colormap(ax,cmap);
        clim(ax,[0.5 4.5]);
        cb = colorbar(ax);
        cb.Ticks = 1:4;
        cb.TickLabels = {'Very low','Low','Medium','High'};
        cb.FontName = fontName;
        cb.FontSize = fontSize;
        cb.TickLength = 0;
        cb.Label.String = 'Cumulative modeled mortality';
        cb.Label.FontName = fontName;
    end

    % Original observed bubbles are shown only on the full baseline map.
    if i == 3
        mort = mortTable{i};
        sizes = nan(height(mort),1);
        ranges = [1 3;4 6;7 10;11 16;17 24;25 33;34 50;51 79];
        for k = 1:size(ranges,1)
            inClass = mort.N_Dead >= ranges(k,1) & mort.N_Dead <= ranges(k,2);
            sizes(inClass) = sizeClass(k);
        end
        keep = isfinite(sizes);
        scatter(ax,mort.X(keep),mort.Y(keep),sizes(keep), ...
            'MarkerEdgeColor','k','MarkerFaceColor','k', ...
            'MarkerFaceAlpha',0.75,'LineWidth',1.0);

        xLeg = 12.5;
        yTop = -23.8;
        dy = 0.58;
        for k = 1:numel(sizeClass)
            yk = yTop-(k-1)*dy;
            scatter(ax,xLeg,yk,sizeClass(k), ...
                'MarkerEdgeColor','k','MarkerFaceColor','k', ...
                'MarkerFaceAlpha',0.75,'LineWidth',1.0);
            if showBubbleLabels
                text(ax,xLeg+0.35,yk,obsLabels{k}, ...
                    'FontSize',20,'FontName',fontName, ...
                    'VerticalAlignment','middle');
            end
        end
    end

    title(ax,titles{i},'FontSize',fontSize,'FontName',fontName, ...
        'FontWeight','normal');
    axis(ax,'equal');
    xlim(ax,xlimMap);
    ylim(ax,ylimMap);
    ax.LineWidth = 1.0;
    ax.FontSize = fontSize;
    ax.FontName = fontName;
    ax.GridAlpha = 0.18;
    ax.TickDir = 'out';
    hold(ax,'off');

    if useInterpolation
        suffix = '_interpolated';
    else
        suffix = '_no_interpolation';
    end
    if showModelBackground
        basename = ['mortality_model_observed_' tags{i} suffix];
    else
        basename = ['observed_mortality_bubbles_' tags{i} suffix];
    end
    exportgraphics(fig,fullfile(outdir,[basename '.eps']),'ContentType','vector');
    exportgraphics(fig,fullfile(outdir,[basename '.png']),'Resolution',300);
    savefig(fig,fullfile(outdir,[basename '.fig']));
    fprintf('Saved %s (.eps, .png, .fig)\n',basename);
end
fprintf('All figures saved to: %s\n',outdir);
