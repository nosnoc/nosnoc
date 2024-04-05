clear all;
[nosnoc_path,~,~] = fileparts(mfilename("fullpath"));
sep = pathsep;
already_installed = contains([sep, path, sep], [sep, nosnoc_path, sep], 'IgnoreCase', ispc);
if ~already_installed
    addpath(nosnoc_path);
    addpath(fullfile(nosnoc_path, 'src')); % TODO: remove once we move everything to +nosnoc
    addpath(fullfile(nosnoc_path, 'external', 'vdx'));
    savepath
end
