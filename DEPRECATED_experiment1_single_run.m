% DEPRECATED: heat-pulse path not pursued (freely moving baseline did not recover loops).
% Kept for reference only; do not use for current experiments.
% Single-worm Experiment 1 run (heat-pulse dataset): fit LOOPER and save results.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % for experiment1_helpers

experiment1_helpers.setup_figs();

% --- Single-worm full E1 settings ---
FIT_BASELINE_ONLY = true;   % fit LOOPER on pre-heat segment only
APPLY_DETREND = false;      % optional linear detrend (fit on pre-heat)

% Load heat worms (per-worm JSONs). Default trace_array is already z-scored.
worms = load_atanas_data('heat');

% Pick the first worm that has a heat event.
eventIdx = find(arrayfun(@(w) isfield(w.events, 'heat_1b'), worms), 1, 'first');
if isempty(eventIdx)
    error('No heat event found in heat split.');
end
worm = worms(eventIdx);

fprintf('Running LOOPER on %s (N=%d, T=%d, heat=%d)\n', ...
    worm.uid, worm.num_neurons, worm.T, worm.events.heat_1b);

% Prepare LOOPER input (neurons x time).
saveData = struct;
rawFull = worm.RawData;
t_on = worm.events.heat_1b;

if FIT_BASELINE_ONLY
    if isempty(worm.ranges_raw)
        error('Expected paper windows (ranges_raw) for Experiment 1.');
    end
    [preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
    experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);
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

saveData.RawData = rawTrain;
saveData.TrialData = ones(1, size(rawTrain, 2)); % single continuous trial
if APPLY_DETREND
    saveData.Detrend = detrendSpec;
else
    saveData.Detrend = struct('enabled', false, 'mode', 'off');
end

params = experiment1_helpers.default_looper_params_atanas();

% Run LOOPER (paper-style settings).
saveData = LOOPER(saveData, true, [], [], [], params);

% Single-worm output bundle.
outDir = fullfile(pwd, 'results', 'experiment1_single');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Save raw final-stream PCA figure from LOOPER (figure 204).
figHandle = findobj('Type', 'figure', 'Number', 204);
if ~isempty(figHandle)
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.fig'));
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.png'));
end

outPath = fullfile(outDir, 'experiment1_single.mat');
save(outPath, 'worm', 'saveData');

fprintf('Saved: %s\\n', outPath);
