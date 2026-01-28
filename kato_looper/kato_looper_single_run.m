% Single-worm LOOPER run on Kato 2015 dataset (immobilized, no-stim).

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % LOOPER.m
addpath(fullfile(rootDir, 'kato_2015'));

% Ensure figures open as docked tabs (not separate windows).
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

% LOOPER parameters: paper-style C. elegans settings (Table 2).
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

worms = load_kato_data();
worm = worms(1);

fprintf('Running LOOPER on Kato worm %s (N=%d, T=%d, dt=%.4fs)\n', ...
    worm.uid, worm.num_neurons, worm.T, worm.dt_sec);

% Prepare LOOPER input (neurons x time).
saveData = struct;
% Detrend per neuron across time (pending OSF clarification).
saveData.RawData = detrend(worm.RawData')';
saveData.TrialData = ones(1, size(worm.RawData, 2));

saveData = LOOPER(saveData, true, [], [], [], params);

% Single-worm output bundle.
outDir = fullfile(rootDir, 'results', 'kato_single');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Save raw final-stream PCA figure from LOOPER (figure 204).
figHandle = findobj('Type', 'figure', 'Number', 204);
if ~isempty(figHandle)
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.fig'));
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.png'));
end

outPath = fullfile(outDir, 'kato_single.mat');
save(outPath, 'worm', 'saveData');

fprintf('Saved: %s\n', outPath);
