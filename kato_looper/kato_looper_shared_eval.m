% Evaluate shared-neuron concatenated Kato LOOPER run.

rootDir = fileparts(fileparts(mfilename('fullpath')));

% Ensure figures open as docked tabs.
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

outDir = fullfile(rootDir, 'results', 'kato_shared');
matPath = fullfile(outDir, 'kato_shared.mat');
if ~exist(matPath, 'file')
    error('No Kato shared results found. Run kato_looper_shared_run first.');
end

S = load(matPath);
if ~isfield(S, 'saveData')
    error('Expected saveData in %s', matPath);
end

diagDir = fullfile(outDir, 'diagnostics');
dt = [];
if isfield(S, 'dt_sec')
    dt = S.dt_sec;
elseif isfield(S.saveData, 'dt_sec')
    dt = S.saveData.dt_sec;
end
diag = run_looper_diagnostics(S.saveData, struct('tag', 'kato_shared', 'outDir', diagDir, 'dt_sec', dt));

meta = struct('uid', 'kato_shared', 'dt_sec', dt);
summary = looper_eval_helpers.base_summary(S.saveData, diag, meta);
writetable(struct2table(summary), fullfile(outDir, 'summary.csv'));

fprintf('Diagnostics saved in: %s\n', diagDir);
