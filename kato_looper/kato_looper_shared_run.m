% Run LOOPER on concatenated Kato worms using the shared neuron set.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % LOOPER.m
addpath(fullfile(rootDir, 'kato_2015'));

% Ensure figures open as docked tabs.
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureVisible', 'on');

params = default_kato_params_();
% We z-score per worm before concatenation; avoid double z-score in LOOPER.
params.PreprocessData.ZScore = false;

outDir = fullfile(rootDir, 'results', 'kato_shared');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
params.CacheDiffusionMapPath = fullfile(outDir, 'diffusion_cache.mat');

% Use identified-label intersection (IDs containing letters).
USE_PAPER_NEURON_LIST = false;

worms = load_kato_data();
if isempty(worms)
    error('No worms found in Kato dataset.');
end
dt_sec = worms(1).dt_sec;

[sharedIds, wormIndices] = shared_neuron_indices_(worms, USE_PAPER_NEURON_LIST);
if isempty(sharedIds)
    error('No shared neuron IDs found across worms.');
end

concatData = [];
trialData = [];
wormOffsets = zeros(1, numel(worms));
wormLengths = zeros(1, numel(worms));
wormUids = strings(1, numel(worms));

for wi = 1:numel(worms)
    worm = worms(wi);
    wormUids(wi) = string(worm.uid);
    idx = wormIndices{wi};

    X = worm.RawData(idx, :); % neurons x time
    % Detrend per neuron across time (pending OSF clarification).
    X = detrend(X')';
    % Z-score per neuron across time (OSF-style).
    X = zscore(X, 0, 2);

    wormOffsets(wi) = size(concatData, 2) + 1;
    wormLengths(wi) = size(X, 2);
    concatData = [concatData, X]; %#ok<AGROW>
    trialData = [trialData, repmat(wi, 1, size(X, 2))]; %#ok<AGROW>
end

saveData = struct;
saveData.RawData = concatData;
saveData.TrialData = trialData;
saveData.dt_sec = dt_sec;

saveData = LOOPER(saveData, true, [], [], [], params);

figHandle = findobj('Type', 'figure', 'Number', 204);
if ~isempty(figHandle)
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.fig'));
    saveas(figHandle, fullfile(outDir, 'final_stream_pca.png'));
end

outPath = fullfile(outDir, 'kato_shared.mat');
save(outPath, 'saveData', 'sharedIds', 'wormOffsets', 'wormLengths', 'wormUids', 'dt_sec');

fprintf('Saved shared-neuron run to: %s\n', outPath);

function params = default_kato_params_()
    params = [];
    params.NearestNeighbors = 8;
    params.UseLocalDimensions = true;
    params.RepopulateDensity = 0.95;
    params.MinimumReturnTime = 50;
    params.DistanceType = 'correlation';
    % Keep MaxCheckTime low to avoid N^2 x MaxCheckTime memory blowups
    % for long concatenated traces (reduceMatrix allocates N x N x T).
    params.MaxCheckTime = 1;
    params.TotalStates = 25;
    params.UseTerminalState = false;
    params.PutativeLoopCounts = [8 7 6 5 4 3 2];
    % Paper-style preprocessing (Table 1).
    params.PreprocessData.ZScore = true;
    params.PreprocessData.Smoothing = 1;
    params.PreprocessData.DelayTime = 10;
    params.PreprocessData.DelayCount = 10;
end

function [sharedIds, wormIndices] = shared_neuron_indices_(worms, usePaperList)
    wormIndices = cell(1, numel(worms));
    allIds = cell(1, numel(worms));

    if usePaperList
        sharedIds = string({ ...
            'AIBL', 'AIBR', 'ALA', 'AVAL', 'AVAR', 'AVBL', 'AVER', 'RID', ...
            'RIML', 'RIMR', 'RMED', 'RMEL', 'RMER', 'VB01', 'VB02' ...
        });
        for wi = 1:numel(worms)
            ids = worms(wi).neuron_ids;
            if isempty(ids)
                error('Worm %d is missing neuron_ids; cannot compute shared set.', wi);
            end
            ids = normalize_ids_(ids);
            [isMember, idx] = ismember(sharedIds, ids);
            if ~all(isMember)
                missing = sharedIds(~isMember);
                error('Worm %d missing paper neuron IDs: %s', wi, strjoin(missing, ', '));
            end
            wormIndices{wi} = idx;
            allIds{wi} = ids;
        end
        return;
    end

    for wi = 1:numel(worms)
        ids = worms(wi).neuron_ids;
        if isempty(ids)
            error('Worm %d is missing neuron_ids; cannot compute shared set.', wi);
        end
        ids = normalize_ids_(ids);
        ids = ids(is_identified_(ids));
        allIds{wi} = ids;
    end

    sharedIds = allIds{1};
    for wi = 2:numel(allIds)
        sharedIds = intersect(sharedIds, allIds{wi}, 'stable');
    end

    for wi = 1:numel(worms)
        [isMember, idx] = ismember(sharedIds, allIds{wi});
        if ~all(isMember)
            error('Shared IDs missing in worm %d after intersection.', wi);
        end
        wormIndices{wi} = idx;
    end
end

    function ids = normalize_ids_(idsIn)
        if isstring(idsIn)
            ids = idsIn;
        elseif ischar(idsIn)
            ids = string(cellstr(idsIn));
        elseif iscell(idsIn)
            ids = string(idsIn);
        else
            ids = string(idsIn);
        end
        ids = ids(:);
    end

    function mask = is_identified_(ids)
        mask = false(size(ids));
        for ii = 1:numel(ids)
            name = char(ids(ii));
            mask(ii) = ~isempty(regexp(name, '[A-Za-z]', 'once'));
        end
    end
