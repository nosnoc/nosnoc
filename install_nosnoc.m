clear all;
[nosnoc_path,~,~] = fileparts(mfilename("fullpath"));
sep = pathsep;
root_already_installed = contains([sep, path, sep], [sep, nosnoc_path, sep], 'IgnoreCase', ispc)
src_already_installed = contains([sep, path, sep], [sep, fullfile(nosnoc_path, 'src'), sep], 'IgnoreCase', ispc)
vdx_already_installed = contains([sep, path, sep], [sep, fullfile(nosnoc_path, 'external', 'vdx'), sep], 'IgnoreCase', ispc)
mpecopt_already_installed = contains([sep, path, sep], [sep, fullfile(nosnoc_path, 'external', 'mpecopt', 'src'), sep], 'IgnoreCase', ispc)

if ~root_already_installed
    addpath(nosnoc_path);
end
if ~src_already_installed
    addpath(fullfile(nosnoc_path, 'src')); % TODO: remove once we move everything to +nosnoc
end
if ~vdx_already_installed
    addpath(fullfile(nosnoc_path, 'external', 'vdx'));
end
if ~mpecopt_already_installed
    addpath(fullfile(nosnoc_path, 'external', 'mpecopt', 'src'));
end
savepath

