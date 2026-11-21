clear; clc; close all;
paths = project_paths();

%% ================= USER SETTINGS =================
dataDir = paths.namibia_data;
baseDir = fullfile(paths.namibia_results, 'full');
outDir  = paths.namibia_figures;

gridCsv  = fullfile(baseDir, "pde_mortality_grid.csv");
mortCsv  = fullfile(dataDir, "Poisoning_GIS.csv");
parksShp = fullfile(dataDir, "Namibia_National_Parks_Boundaries_2001.shp");

overlayObservedPoisoning = true;
overlayParks = true;

xlimMap = [11.99489, 21.0];
ylimMap = [-29.0, -16.9059];

cellSize = 20;
smoothSigma = 0.5;

%% ================= LOAD GRID =================
G = readtable(gridCsv);

lon = G.x;
lat = G.y;
food = G.food_model;
toxin = G.toxin_model;

valid = isfinite(lon) & isfinite(lat) & isfinite(food) & isfinite(toxin);

lon = lon(valid);
lat = lat(valid);
food = food(valid);
toxin = toxin(valid);

%% ================= OPTIONAL MORTALITY =================
M = table();        % important: initialize M
hasMort = false;

if overlayObservedPoisoning && isfile(mortCsv)
    M = readtable(mortCsv);

    if all(ismember(["X","Y","N_Dead","Death_Reason"], string(M.Properties.VariableNames)))
        reason = lower(strtrim(string(M.Death_Reason)));

        keep = (reason == "poisoning" | reason == "probably_poisoning") ...
            & isfinite(M.X) & isfinite(M.Y) & M.N_Dead > 0;

        M = M(keep,:);
        hasMort = height(M) > 0;
    end
else
    warning("Mortality file not found: %s", mortCsv);
end

%% ================= SMOOTH =================
food_smooth  = imgaussfilt(food, smoothSigma);
toxin_smooth = imgaussfilt(toxin, smoothSigma);

%% ================= NORMALIZE =================
food_smooth = food_smooth - min(food_smooth);
food_smooth = food_smooth / max(food_smooth);

toxin_smooth = toxin_smooth - min(toxin_smooth);
toxin_smooth = toxin_smooth / max(toxin_smooth);

%% =================================================
%% FIGURE 1 — FOOD
%% =================================================
gammaFood = 0.45;
food_display = gammaColorStretch(food_smooth, gammaFood);

fig1 = figure("Color","w");
ax1 = axes(fig1);

scatter(ax1, lon, lat, cellSize, food_display, "s", "filled", ...
    "MarkerEdgeColor","none");

hold(ax1,"on");

decorate(ax1, xlimMap, ylimMap);

colormap(ax1, cool);
clim(ax1,[0 1]);

cb1 = colorbar(ax1);
cb1.FontSize = 18;

tickVals = [0 0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8 1.0];
cb1.Ticks = gammaColorStretch(tickVals, gammaFood);
cb1.TickLabels = ["0","0.01","0.02","0.05","0.1","0.2","0.4","0.6","0.8","1"];

xlabel(ax1,"Longitude");
ylabel(ax1,"Latitude");

overlay(ax1, overlayParks, parksShp, M, hasMort);

%% =================================================
%% FIGURE 2 — TOXIN
%% =================================================
fig2 = figure("Color","w");
ax2 = axes(fig2);

scatter(ax2, lon, lat, cellSize, toxin_smooth, "s", "filled", ...
    "MarkerEdgeColor","none");

hold(ax2,"on");

decorate(ax2, xlimMap, ylimMap);

colormap(ax2, hot(256));
clim(ax2,[0 1]);

cb2 = colorbar(ax2);
cb2.FontSize = 18;

xlabel(ax2,"Longitude");
ylabel(ax2,"Latitude");

 overlay(ax2, overlayParks, parksShp, M, hasMort);

%% ================= SAVE =================
if ~isfolder(outDir)
    mkdir(outDir);
end

exportgraphics(fig1, fullfile(outDir,"food_map.png"), "Resolution", 300);
exportgraphics(fig2, fullfile(outDir,"toxin_map.png"), "Resolution", 300);
exportgraphics(fig1, fullfile(outDir,"food_map.eps"), "ContentType", "vector");
exportgraphics(fig2, fullfile(outDir,"toxin_map.eps"), "ContentType", "vector");
savefig(fig1, fullfile(outDir,"food_map.fig"));
savefig(fig2, fullfile(outDir,"toxin_map.fig"));

fprintf("Saved figures to:\n%s\n", outDir);

%% =================================================
%% LOCAL FUNCTIONS
%% =================================================
function decorate(ax, xlimMap, ylimMap)
    axis(ax,"equal");
    axis(ax,"tight");

    xlim(ax,xlimMap);
    ylim(ax,ylimMap);

    box(ax,"on");
    grid(ax,"on");

    ax.LineWidth = 1.0;
    ax.FontSize = 18;
    ax.GridAlpha = 0.18;
    ax.TickDir = "out";
end

function overlay(ax, overlayParks, parksShp, M, hasMort)

    % --- parks ---
    if overlayParks && isfile(parksShp) && exist("shaperead","file") == 2
        try
            S = shaperead(parksShp);
            for k = 1:numel(S)
                plot(ax, S(k).X, S(k).Y, ...
                    "Color",[0.2 0.2 0.2], ...
                    "LineWidth",0.7);
            end
        catch
        end
    end

    % --- observed mortality ---
    if hasMort
        scatter(ax, M.X, M.Y, sqrt(M.N_Dead)*25, ...
            "o", "MarkerEdgeColor","k", ...
            "MarkerFaceColor","none", "LineWidth",0.9);
    end
end

function y = gammaColorStretch(x, gammaValue)
    x = max(0, min(1, x));
    y = x .^ gammaValue;
end


