% Evaluate Kato LOOPER run using built-in diagnostics.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % run_looper_diagnostics
addpath(fullfile(rootDir, 'kato_2015'));

% Ensure figures open as docked tabs (not separate windows).
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

outDir = fullfile(rootDir, 'results', 'kato_single');
matPath = fullfile(outDir, 'kato_single.mat');
if ~exist(matPath, 'file')
    error('No Kato single-worm results found. Run kato_looper_single_run first.');
end

S = load(matPath);
if ~isfield(S, 'saveData')
    error('Expected saveData in %s', matPath);
end

saveData = S.saveData;
tag = 'kato';
if isfield(S, 'worm') && isfield(S.worm, 'uid')
    tag = sprintf('kato_worm_%s', S.worm.uid);
end

diagDir = fullfile(fileparts(matPath), 'diagnostics');
dt = [];
if isfield(S, 'worm') && isfield(S.worm, 'dt_sec')
    dt = S.worm.dt_sec;
end
diag = run_looper_diagnostics(saveData, struct('tag', tag, 'outDir', diagDir, 'dt_sec', dt));

meta = struct('uid', tag, 'dt_sec', dt);
if isfield(S, 'worm')
    if isfield(S.worm, 'uid'); meta.uid = sprintf('kato_worm_%s', S.worm.uid); end
    if isfield(S.worm, 'num_neurons'); meta.n_neurons = S.worm.num_neurons; end
    if isfield(S.worm, 'T'); meta.T = S.worm.T; end
end
summary = looper_eval_helpers.base_summary(saveData, diag, meta);
writetable(struct2table(summary), fullfile(outDir, 'summary.csv'));

figPath = fullfile(fileparts(matPath), 'final_stream_pca.fig');
if exist(figPath, 'file')
    openfig(figPath, 'new', 'visible');
end

fprintf('Diagnostics saved in: %s\n', diagDir);
