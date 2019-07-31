function GEOREPRStartUp()

% This function loads FVTool and the folders of original submodels

p = mfilename('fullpath');
file_name = mfilename;
current_path = p(1:end-1-length(file_name));
addpath([current_path '/FVTool']);
addpath([current_path '/Utilities']);
addpath([current_path '/Reactions']);
addpath([current_path '/Sourceterms']);

FVToolStartUp()

end
