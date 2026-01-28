% Run LOOPER across all Atanas baseline worms and save per-worm results.

addpath(fullfile(pwd, 'atanas-data'));
addpath(pwd); % experiment1_helpers

experiment1_helpers.setup_figs();

% --- Settings (match experiment1_single_baseline_run) ---
FIT_BASELINE_ONLY = true; % fit LOOPER on first half only
APPLY_DETREND = false;    % optional linear detrend (fit on pre)

params = experiment1_helpers.default_looper_params_atanas();

worms = load_atanas_data('baseline');

outDir = fullfile(pwd, 'results', 'experiment1_all_baseline');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for i = 1:numel(worms)
    worm = worms(i);

    if isempty(worm.ranges_raw)
        fprintf('Skipping %s (no ranges_raw).\n', worm.uid);
        continue;
    end

    [preIdx, postIdx] = experiment1_helpers.normalize_ranges(worm.ranges_raw);
    t_on = postIdx(1); %#ok<NASGU>
    experiment1_helpers.warn_ranges_contiguous(preIdx, postIdx, worm.uid);

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

    fprintf('Running LOOPER on baseline %s (N=%d, T=%d)\n', ...
        worm.uid, worm.num_neurons, worm.T);

    safeUid = regexprep(worm.uid, '[^A-Za-z0-9_-]', '_');
    cacheDir = fullfile(outDir, 'cache');
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end
    params.CacheDiffusionMapPath = fullfile(cacheDir, sprintf('worm_%s_diffusion_cache.mat', safeUid));

    saveData = LOOPER(saveData, true, [], [], [], params);

    % Save raw final-stream PCA figure from LOOPER (figure 204).
    figHandle = findobj('Type', 'figure', 'Number', 204);
    if ~isempty(figHandle)
        saveas(figHandle, fullfile(outDir, sprintf('worm_%s_final_stream_pca.fig', safeUid)));
        saveas(figHandle, fullfile(outDir, sprintf('worm_%s_final_stream_pca.png', safeUid)));
    end

    outPath = fullfile(outDir, sprintf('worm_%s.mat', safeUid));
    save(outPath, 'worm', 'saveData');

    fprintf('Saved: %s\n', outPath);
end
