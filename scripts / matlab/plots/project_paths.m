function paths = project_paths()
%PROJECT_PATHS Resolve input and output directories for MATLAB plots.
% Run any plot from any working directory after adding this folder to the
% MATLAB path. Optional environment variables override the defaults.

root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
paths.root = root;
paths.namibia_data = directory('NAMIBIA_DATA_DIR', fullfile(root, 'data', 'namibia'));
paths.namibia_results = directory('NAMIBIA_RESULTS_DIR', fullfile(root, 'results', 'namibia'));
paths.namibia_figures = directory('NAMIBIA_FIGURES_DIR', fullfile(root, 'figures', 'namibia'));
paths.agent_results = directory('AGENT_RESULTS_DIR', fullfile(root, 'results', 'agent'));
paths.agent_figures = directory('AGENT_FIGURES_DIR', fullfile(root, 'figures', 'agent'));
end

function value = directory(variable, fallback)
value = getenv(variable);
if isempty(value)
    value = fallback;
end
end
