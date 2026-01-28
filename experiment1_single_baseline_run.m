% Single-worm baseline run: fit LOOPER on a no-heat worm and save results.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % for experiment1_helpers

experiment1_helpers.setup_figs();

% --- Settings (same as experiment1_single_run) ---
FIT_BASELINE_ONLY = true;   % fit LOOPER on first half only
APPLY_DETREND = false;       % optional linear detrend (fit on pre)

params = experiment1_helpers.default_looper_params_atanas();

% Load baseline worms (per-worm JSONs). Default trace_array is already z-scored.
worms = load_atanas_data('baseline');
worm = worms(1);

% Use the split between ranges as a pseudo-event.
if isempty(worm.ranges_raw)
    error('Expected paper windows (ranges_raw) for Experiment 1 baseline.');
end
[preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
t_on = postIdx(1); % boundary between halves
experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);

fprintf('Running LOOPER on baseline %s (N=%d, T=%d, split=%d)\n', ...
    worm.uid, worm.num_neurons, worm.T, t_on);

% Fit LOOPER on pre segment only.
rawFull = worm.RawData;
if FIT_BASELINE_ONLY
    if APPLY_DETREND
        detrendSpec = experiment1_helpers.detrend_fit(rawFull, preIdx);
        rawFull = experiment1_helpers.detrend_apply(rawFull, detrendSpec);
    end
    rawTrain = rawFull(:, preIdx);
else
    if APPLY_DETREND
        detrendSpec = experiment1_helpers.detrend_fit(rawFull, []);
        rawFull = experiment1_helpers.detrend_apply(rawFull, detrendSpec);
    end
    rawTrain = rawFull;
end

saveData = struct;
saveData.RawData = rawTrain;
saveData.TrialData = ones(1, size(rawTrain, 2));
if APPLY_DETREND
    saveData.Detrend = detrendSpec;
else
    saveData.Detrend = struct('enabled', false, 'mode', 'off');
end

saveData = LOOPER(saveData, true, [], [], [], params);

% Save bundle.
outDir = fullfile(pwd, 'results', 'experiment1_single_baseline');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Save raw final-stream PCA figure from LOOPER (figure 204).
figHandle = findobj('Type', 'figure', 'Number', 204);
if ~isempty(figHandle)
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.fig'));
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.png'));
end

outPath = fullfile(outDir, 'experiment1_single_baseline.mat');
save(outPath, 'worm', 'saveData');

fprintf('Saved: %s\n', outPath);
