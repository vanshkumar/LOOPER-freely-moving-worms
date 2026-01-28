% Run LOOPER diagnostics on the OSF wormlooper2.mat checkpoint.
%
% Notes:
% - This is for inspecting the published OSF checkpoint, not for new fits.
% - `run_looper_diagnostics` has guards for variable trial lengths and will
%   skip validateData when the OSF assumptions do not hold.
% - Output is written to results/osf_wormlooper2_diagnostics.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % run_looper_diagnostics

% Ensure figures open as docked tabs.
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

matPath = fullfile(rootDir, 'LOOPER_OSF_hosted_2022', 'wormlooper2.mat');
if ~exist(matPath, 'file')
    error('Missing wormlooper2.mat at %s', matPath);
end

outDir = fullfile(rootDir, 'results', 'osf_wormlooper2_diagnostics');
run_looper_diagnostics(matPath, struct('tag', 'wormlooper2', 'outDir', outDir));

fprintf('Diagnostics saved in: %s\n', outDir);
