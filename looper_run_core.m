function looper_run_core(dataset, mode, scope)
%LOOPER_RUN_CORE Shared runner for Atanas + Kato LOOPER runs.
%   dataset: "atanas" | "kato"
%   mode:    "stationarity" | "fidelity"
%   scope:   "single" | "all"

    if nargin < 3
        error('looper_run_core requires dataset, mode, and scope.');
    end
    dataset = lower(string(dataset));
    mode = lower(string(mode));
    scope = lower(string(scope));
    if dataset ~= "atanas" && dataset ~= "kato"
        error('dataset must be "atanas" or "kato".');
    end
    if mode ~= "stationarity" && mode ~= "fidelity"
        error('mode must be "stationarity" or "fidelity".');
    end
    if scope ~= "single" && scope ~= "all"
        error('scope must be "single" or "all".');
    end

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    if dataset == "atanas"
        addpath(fullfile(rootDir, 'atanas-data'));
    else
        addpath(fullfile(rootDir, 'kato_2015'));
        addpath(fullfile(rootDir, 'kato_looper'));
    end

    set(0, 'DefaultFigureVisible', 'off');

    fitStationarity = mode == "stationarity";
    if dataset == "atanas"
        params = atanas_helpers.default_looper_params_atanas();
        APPLY_DETREND = false;
        worms = load_atanas_data('baseline');
    else
        params = kato_helpers.default_kato_params();
        APPLY_DETREND = true;
        worms = load_kato_data();
    end

    if isempty(worms)
        error('No worms found for dataset %s.', dataset);
    end

    outDir = results_dir_(rootDir, dataset, scope, mode);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    if scope == "single"
        worm = worms(1);
        fprintf('Running LOOPER on %s worm %s (%s) (N=%d, T=%d, dt=%.4fs)\n', ...
            dataset, worm.uid, mode, worm.num_neurons, worm.T, worm.dt_sec);

        saveData = struct;
        [rawTrain, detrendSpec, split] = prepare_training_(dataset, worm, fitStationarity, APPLY_DETREND);
        saveData.RawData = rawTrain;
        saveData.TrialData = ones(1, size(rawTrain, 2)); % single continuous trial
        saveData.Detrend = detrendSpec;
        saveData.TrainSplit = split;

        saveData = LOOPER(saveData, true, [], [], [], params);
        % LOOPER overwrites saveData; reattach run metadata needed by eval.
        saveData.Detrend = detrendSpec;
        saveData.TrainSplit = split;

        outPath = fullfile(outDir, sprintf('%s_%s_%s.mat', dataset, scope, mode));
        save(outPath, 'worm', 'saveData');
        fprintf('Saved: %s\n', outPath);
        return;
    end

    cacheDir = fullfile(outDir, 'cache');
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end

    manifest = struct;
    manifest.params = params;
    manifest.worm_uids = strings(1, numel(worms));
    manifest.mat_paths = strings(1, numel(worms));

    for wi = 1:numel(worms)
        worm = worms(wi);
        fprintf('Running LOOPER on %s worm %s (%s) (N=%d, T=%d, dt=%.4fs)\n', ...
            dataset, worm.uid, mode, worm.num_neurons, worm.T, worm.dt_sec);

        saveData = struct;
        [rawTrain, detrendSpec, split] = prepare_training_(dataset, worm, fitStationarity, APPLY_DETREND);
        saveData.RawData = rawTrain;
        saveData.TrialData = ones(1, size(rawTrain, 2)); % single continuous trial
        saveData.Detrend = detrendSpec;
        saveData.TrainSplit = split;

        safeUid = regexprep(worm.uid, '[^A-Za-z0-9_-]', '_');
        params.CacheDiffusionMapPath = fullfile(cacheDir, sprintf('worm_%s_diffusion_cache.mat', safeUid));

        saveData = LOOPER(saveData, true, [], [], [], params);
        % LOOPER overwrites saveData; reattach run metadata needed by eval.
        saveData.Detrend = detrendSpec;
        saveData.TrainSplit = split;

        outPath = fullfile(outDir, sprintf('worm_%s.mat', safeUid));
        save(outPath, 'worm', 'saveData');

        manifest.worm_uids(wi) = string(worm.uid);
        manifest.mat_paths(wi) = string(outPath);
    end

    save(fullfile(outDir, 'manifest.mat'), 'manifest');
    fprintf('Saved all-worm results in: %s\n', outDir);
end

function outDir = results_dir_(rootDir, dataset, scope, mode)
    outDir = fullfile(rootDir, 'results', sprintf('%s_%s', dataset, scope), mode);
end

function [rawTrain, detrendSpec, split] = prepare_training_(dataset, worm, fitStationarity, applyDetrend)
    if dataset == "atanas"
        [rawTrain, detrendSpec, split] = atanas_helpers.prepare_training(worm, fitStationarity, applyDetrend);
        return;
    end
    [rawTrain, detrendSpec, split] = kato_helpers.prepare_training(worm, fitStationarity, applyDetrend);
end
