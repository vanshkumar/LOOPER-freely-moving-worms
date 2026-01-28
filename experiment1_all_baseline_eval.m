% Evaluate Experiment 1 baseline runs across all worms.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % experiment1_helpers, run_looper_diagnostics

% Ensure figures open as docked tabs (not separate windows).
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

RUN_DIAGNOSTICS = true;

outDir = fullfile(pwd, 'results', 'experiment1_all_baseline');
if ~exist(outDir, 'dir')
    error('No results folder found: %s. Run experiment1_all_baseline_run first.', outDir);
end

files = dir(fullfile(outDir, 'worm_*.mat'));
if isempty(files)
    error('No baseline results found in %s. Run experiment1_all_baseline_run first.', outDir);
end

plotDir = fullfile(outDir, 'plots');
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

projDir = fullfile(outDir, 'projections');
if ~exist(projDir, 'dir')
    mkdir(projDir);
end

diagRoot = fullfile(outDir, 'diagnostics');
if RUN_DIAGNOSTICS && ~exist(diagRoot, 'dir')
    mkdir(diagRoot);
end

summary = struct('uid', {}, 'n_neurons', {}, 'T', {}, 'dt_sec', {}, ...
    'split_idx', {}, 'unique_loops', {}, 'loop_switches', {}, 'segments', {}, ...
    'mean_pre', {}, 'mean_post', {}, 'post_slope', {}, 'd_peak', {}, ...
    'phase_frac_small', {}, 'phase_var', {}, 'cycles_per_min', {}, ...
    'median_d', {}, 'mean_segment_len', {});

for i = 1:numel(files)
    S = load(fullfile(files(i).folder, files(i).name));
    if ~isfield(S, 'worm') || ~isfield(S, 'saveData')
        fprintf('Skipping %s (missing worm/saveData).\n', files(i).name);
        continue;
    end

    worm = S.worm;
    saveData = S.saveData;

    if isempty(worm.ranges_raw)
        fprintf('Skipping %s (no ranges_raw).\n', worm.uid);
        continue;
    end

    [preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
    t_on = postIdx(1);
    experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);

    rawFull = worm.RawData;
    if isfield(saveData, 'Detrend') && isfield(saveData.Detrend, 'enabled') && saveData.Detrend.enabled
        rawFull = experiment1_helpers.detrend_apply(rawFull, saveData.Detrend);
    end
    [alpha, theta, stateIdx, d] = experiment1_helpers.project_to_model(rawFull, saveData); %#ok<ASGLU>

    

    dt = worm.dt_sec;
    if isempty(dt) || ~isfinite(dt)
        dt = 1;
    end

    out = experiment1_helpers.compute_delta_d(d, dt, saveData, t_on, preIdx, postIdx);
    relTime = out.relTime;
    delta_d = out.delta_d;
    tOnProc = out.tOnProc;

    safeUid = regexprep(worm.uid, '[^A-Za-z0-9_-]', '_');

    % Baseline recovery plot.
    figRecovery = figure('Name', sprintf('Baseline recovery (%s)', worm.uid), ...
        'NumberTitle', 'off');
    plot(relTime, delta_d, 'LineWidth', 1.5);
    xline(0, '--k');
    xlabel('Seconds from split');
    ylabel('\Delta d (baseline-normalized)');
    title(sprintf('Baseline recovery (worm %s)', worm.uid));
    saveas(figRecovery, fullfile(plotDir, sprintf('worm_%s_baseline_recovery.png', safeUid)));

    % Loop assignment over time (embedded timebase).
    figLoops = figure('Name', sprintf('Loop assignment (%s)', worm.uid), ...
        'NumberTitle', 'off');
    plot(alpha, 'LineWidth', 1.2);
    xline(tOnProc, '--k');
    xlabel('Time (embedded frames)');
    ylabel('Loop assignment');
    title(sprintf('Loop assignment over time (worm %s)', worm.uid));
    saveas(figLoops, fullfile(plotDir, sprintf('worm_%s_loop_assignment.png', safeUid)));

    save(fullfile(projDir, sprintf('worm_%s_projection.mat', safeUid)), ...
        'alpha', 'theta', 'stateIdx', 'd');

    % Summary metrics.
    nBins = max(saveData.BestLoopAssignments(:,2));
    if isempty(nBins) || ~isfinite(nBins)
        nBins = max(theta);
    end
    summaryRow = experiment1_helpers.summary_metrics( ...
        alpha, theta, d, relTime, delta_d, worm, t_on, dt, nBins);

    summary(end+1) = summaryRow; %#ok<SAGROW>

    if RUN_DIAGNOSTICS
        diagDir = fullfile(diagRoot, safeUid);
        run_looper_diagnostics(saveData, struct('tag', safeUid, 'outDir', diagDir, ...
            'dt_sec', dt, 'preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));
    end
end

if ~isempty(summary)
    T = struct2table(summary);
    writetable(T, fullfile(outDir, 'summary.csv'));
end

fprintf('Baseline eval complete. Summary saved to %s\n', fullfile(outDir, 'summary.csv'));
