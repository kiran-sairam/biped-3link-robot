% set_path.m

project_root = fileparts(mfilename('fullpath'));

addpath(project_root);
addpath(fullfile(project_root, 'util'));
addpath(fullfile(project_root, 'autogen'));

disp('Project paths added successfully.');