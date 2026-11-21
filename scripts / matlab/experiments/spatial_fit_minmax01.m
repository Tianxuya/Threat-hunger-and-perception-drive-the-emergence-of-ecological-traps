clear; clc; close all;

% The baseline PDE run writes mortality_event_comparison.csv here.
script_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(fileparts(script_dir)));
namibia_results_dir = getenv('NAMIBIA_RESULTS_DIR');
if isempty(namibia_results_dir), namibia_results_dir = fullfile(repo_dir, 'results', 'namibia'); end
outdir = fullfile(namibia_results_dir, 'full');
figures_dir = getenv('NAMIBIA_FIGURES_DIR');
if isempty(figures_dir), figures_dir = fullfile(repo_dir, 'figures', 'namibia'); end
if ~isfolder(figures_dir), mkdir(figures_dir); end
if ~isfile(fullfile(outdir, 'mortality_event_comparison.csv'))
    error('Missing %s. Run the baseline Namibia PDE model first.', ...
        fullfile(outdir, 'mortality_event_comparison.csv'));
end

mort = readtable( ...
    fullfile(outdir, "mortality_event_comparison.csv"));
required_columns = {'modeled_total_mortality_mass_within_radius_minmax01', ...
    'observed_deaths_within_radius_minmax01'};
if ~all(ismember(required_columns, mort.Properties.VariableNames))
    error('The mortality CSV needs columns: %s.', strjoin(required_columns, ', '));
end

disp(mort.Properties.VariableNames')

% Use min-max normalized columns
x = mort.modeled_total_mortality_mass_within_radius_minmax01;
y = mort.observed_deaths_within_radius_minmax01;

valid = isfinite(x) & isfinite(y);

x = x(valid);
y = y(valid);
if numel(x) < 2
    error('At least two finite mortality pairs are required.');
end

pearson_r = corr(x, y, "Type", "Pearson");
spearman_r = corr(x, y, "Type", "Spearman");

resid = x - y;

mae_val  = mean(abs(resid));
mse_val  = mean(resid.^2);
rmse_val = sqrt(mse_val);

figure( ...
    "Color", "w", ...
    "Position", [100 100 900 650]);

scatter( ...
    x, y, 55, "filled", ...
    "MarkerFaceAlpha", 0.65, ...
    "MarkerEdgeColor", [0.1 0.3 0.5]);

xlabel( ...
    "Modeled total mortality within 50 km, scaled 0-1", ...
    "FontSize", 14);

ylabel( ...
    "Observed vulture deaths within 50 km, scaled 0-1", ...
    "FontSize", 14);

title( ...
    "Min-max normalized spatial fit of mortality", ...
    "FontSize", 16);

txt = sprintf( ...
    "Pearson = %.3f\nSpearman = %.3f\nMAE = %.4f\nMSE = %.4f\nRMSE = %.4f", ...
    pearson_r, ...
    spearman_r, ...
    mae_val, ...
    mse_val, ...
    rmse_val);

text( ...
    0.05, 0.95, txt, ...
    "Units", "normalized", ...
    "VerticalAlignment", "top", ...
    "FontSize", 14);

axis([-0.05 1.05 -0.05 1.05]);

box on;
grid off;

exportgraphics( ...
    gcf, ...
    fullfile(figures_dir, "matlab_minmax01_spatial_fit.png"), ...
    "Resolution", 300);
