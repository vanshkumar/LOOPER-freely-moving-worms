% Run LOOPER on all 5 Kato worms (independent fits).

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % LOOPER.m
addpath(fullfile(rootDir, 'kato_2015'));

% Ensure figures open as docked tabs.
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

params = default_kato_params_();

worms = load_kato_data();
outDir = fullfile(rootDir, 'results', 'kato_all');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

manifest = struct;
manifest.params = params;
manifest.worm_uids = strings(1, numel(worms));
manifest.mat_paths = strings(1, numel(worms));

for wi = 1:numel(worms)
    worm = worms(wi);
    fprintf('Running LOOPER on Kato worm %s (N=%d, T=%d, dt=%.4fs)\n', ...
        worm.uid, worm.num_neurons, worm.T, worm.dt_sec);

    saveData = struct;
    % Detrend per neuron across time (pending OSF clarification).
    saveData.RawData = detrend(worm.RawData')';
    saveData.TrialData = ones(1, size(worm.RawData, 2));

    saveData = LOOPER(saveData, true, [], [], [], params);

    figHandle = findobj('Type', 'figure', 'Number', 204);
    if ~isempty(figHandle)
        saveas(figHandle, fullfile(outDir, sprintf('worm_%s_final_stream_pca.fig', worm.uid)));
        saveas(figHandle, fullfile(outDir, sprintf('worm_%s_final_stream_pca.png', worm.uid)));
    end

    outPath = fullfile(outDir, sprintf('worm_%s.mat', worm.uid));
    save(outPath, 'worm', 'saveData');

    manifest.worm_uids(wi) = string(worm.uid);
    manifest.mat_paths(wi) = string(outPath);
end

save(fullfile(outDir, 'manifest.mat'), 'manifest');
fprintf('Saved all-worm results in: %s\n', outDir);

function params = default_kato_params_()
    params = [];
    params.NearestNeighbors = 8;
    params.UseLocalDimensions = true;
    params.RepopulateDensity = 0.95;
    params.MinimumReturnTime = 10;
    params.DistanceType = 'correlation';
    params.MaxCheckTime = 10;
    params.TotalStates = 25;
    params.UseTerminalState = false;
    params.PutativeLoopCounts = [8 7 6 5 4 3 2];
    % Paper-style preprocessing (Table 1). Kato traces are not pre-zscored.
    params.PreprocessData.ZScore = true;
    params.PreprocessData.Smoothing = 1;
    params.PreprocessData.DelayTime = 10;
    params.PreprocessData.DelayCount = 10;
end
