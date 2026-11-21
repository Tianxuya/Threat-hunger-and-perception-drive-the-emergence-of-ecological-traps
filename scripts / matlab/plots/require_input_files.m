function require_input_files(folder, filenames)
%REQUIRE_INPUT_FILES Report missing precomputed simulation inputs together.
missing = {};
for i = 1:numel(filenames)
    filename = fullfile(folder, filenames{i});
    if ~isfile(filename)
        missing{end+1} = filename; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('Missing required input files:\n%s', strjoin(missing, newline));
end
end
