% Evaluate all-worm Kato LOOPER runs using built-in diagnostics.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % run_looper_diagnostics
addpath(fullfile(rootDir, 'kato_2015'));

% Ensure figures open as docked tabs.
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

outDir = fullfile(rootDir, 'results', 'kato_all');
files = dir(fullfile(outDir, 'worm_*.mat'));
if isempty(files)
    error('No Kato all-worm results found. Run kato_looper_all_run first.');
end

diagRoot = fullfile(outDir, 'diagnostics');
if ~exist(diagRoot, 'dir')
    mkdir(diagRoot);
end

summaries = struct('uid', {}, 'n_neurons', {}, 'T', {}, 'dt_sec', {}, ...
    'unique_loops', {}, 'loop_switches', {}, 'segments', {}, ...
    'validation_score_mean', {}, 'validation_score_std', {}, ...
    'phase_frac_small', {}, 'phase_var', {}, 'cycles_per_min', {}, ...
    'median_d', {}, 'mean_segment_len', {});

for i = 1:numel(files)
    matPath = fullfile(files(i).folder, files(i).name);
    S = load(matPath);
    if ~isfield(S, 'saveData')
        warning('Skipping %s (missing saveData).', matPath);
        continue;
    end

    tag = 'kato_worm';
    if isfield(S, 'worm') && isfield(S.worm, 'uid')
        tag = sprintf('kato_worm_%s', S.worm.uid);
    end

    diagDir = fullfile(diagRoot, tag);
    dt = [];
    if isfield(S, 'worm') && isfield(S.worm, 'dt_sec')
        dt = S.worm.dt_sec;
    end
    diag = run_looper_diagnostics(S.saveData, struct('tag', tag, 'outDir', diagDir, 'dt_sec', dt));

    meta = struct('uid', tag, 'dt_sec', dt);
    if isfield(S, 'worm')
        if isfield(S.worm, 'uid'); meta.uid = sprintf('kato_worm_%s', S.worm.uid); end
        if isfield(S.worm, 'num_neurons'); meta.n_neurons = S.worm.num_neurons; end
        if isfield(S.worm, 'T'); meta.T = S.worm.T; end
    end
    summaries(end+1) = looper_eval_helpers.base_summary(S.saveData, diag, meta); %#ok<SAGROW>
end

if ~isempty(summaries)
    writetable(struct2table(summaries), fullfile(outDir, 'summary.csv'));
end

fprintf('Diagnostics saved in: %s\n', diagRoot);
