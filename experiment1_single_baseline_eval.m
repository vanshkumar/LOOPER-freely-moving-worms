% Evaluate Experiment 1 single baseline run: recovery plot + LOOPER diagnostics.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % experiment1_helpers, run_looper_diagnostics

% Ensure figures open as docked tabs (not separate windows).
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

% --- Settings (match run) ---

outDir = fullfile(pwd, 'results', 'experiment1_single_baseline');
matPath = fullfile(outDir, 'experiment1_single_baseline.mat');
if ~exist(matPath, 'file')
    error('No Experiment1 single baseline results found. Run experiment1_single_baseline_run first.');
end

S = load(matPath);
if ~isfield(S, 'saveData') || ~isfield(S, 'worm')
    error('Expected worm and saveData in %s', matPath);
end

worm = S.worm;
saveData = S.saveData;

% Use the split between ranges as a pseudo-event.
if isempty(worm.ranges_raw)
    error('Expected paper windows (ranges_raw) for Experiment 1 baseline.');
end
[preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
t_on = postIdx(1);
experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);

% Project full trace onto learned scaffold.
rawFull = worm.RawData;
if isfield(saveData, 'Detrend') && isfield(saveData.Detrend, 'enabled') && saveData.Detrend.enabled
    rawFull = experiment1_helpers.detrend_apply(rawFull, saveData.Detrend);
end
[alpha, theta, stateIdx, d] = experiment1_helpers.project_to_model(rawFull, saveData);

% Plot delta d around the split (delay-embedded timebase).
dt = worm.dt_sec;
out = experiment1_helpers.compute_delta_d(d, dt, saveData, t_on, preIdx, postIdx);
relTime = out.relTime;
delta_d = out.delta_d;

outDir = fullfile(pwd, 'results', 'experiment1_single_baseline');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figure(11); clf;
plot(relTime, delta_d, 'LineWidth', 1.5);
xlabel('Seconds from split');
ylabel('\Delta d (baseline-normalized)');
title(sprintf('Baseline recovery (worm %s)', worm.uid));
saveas(gcf, fullfile(outDir, 'baseline_recovery.png'));

save(fullfile(outDir, 'baseline_projection.mat'), ...
    'alpha', 'theta', 'stateIdx', 'd');

% Summary CSV (single-row, all-baseline-compatible fields).
nBins = max(saveData.BestLoopAssignments(:,2));
if isempty(nBins) || ~isfinite(nBins)
    nBins = max(theta);
end
summary = experiment1_helpers.summary_metrics( ...
    alpha, theta, d, relTime, delta_d, worm, t_on, dt, nBins);
summaryPath = fullfile(outDir, 'summary.csv');
writetable(struct2table(summary), summaryPath);

% --- Diagnostics (LOOPER scripts) ---
diagDir = fullfile(outDir, 'diagnostics');
diagTag = sprintf('atanas_baseline_%s', worm.uid);
run_looper_diagnostics(saveData, struct('tag', diagTag, 'outDir', diagDir, ...
    'dt_sec', dt, 'preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));

figPath = fullfile(outDir, 'final_stream_pca.fig');
if exist(figPath, 'file')
    openfig(figPath, 'new', 'visible');
end

fprintf('Diagnostics saved in: %s\n', diagDir);
