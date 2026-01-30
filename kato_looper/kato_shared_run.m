% Run LOOPER on concatenated Kato worms using the shared neuron set.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir); % LOOPER.m
addpath(fullfile(rootDir, 'kato_2015'));
addpath(fullfile(rootDir, 'kato_looper'));

set(0, 'DefaultFigureVisible', 'off');

params = kato_helpers.default_kato_params();
% We z-score per worm before concatenation; avoid double z-score in LOOPER.
params.PreprocessData.ZScore = false;
% OSF worm checkpoint uses MinimumReturnTime = 50 (vs paper default 10).
% We keep that here to match the shared-worm run behavior.
params.MinimumReturnTime = 50;
% Keep MaxCheckTime low to avoid N^2 x MaxCheckTime memory blowups
% for long concatenated traces (reduceMatrix allocates N x N x T).
params.MaxCheckTime = 1;

outDir = fullfile(rootDir, 'results', 'kato_shared');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
params.CacheDiffusionMapPath = fullfile(outDir, 'diffusion_cache.mat');

% Use identified-label intersection (IDs containing letters).

worms = load_kato_data();
if isempty(worms)
    error('No worms found in Kato dataset.');
end
dt_sec = worms(1).dt_sec;

[sharedIds, wormIndices] = shared_neuron_indices_(worms);
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

outPath = fullfile(outDir, 'kato_shared.mat');
save(outPath, 'saveData', 'sharedIds', 'wormOffsets', 'wormLengths', 'wormUids', 'dt_sec');

fprintf('Saved shared-neuron run to: %s\n', outPath);

function [sharedIds, wormIndices] = shared_neuron_indices_(worms)
    wormIndices = cell(1, numel(worms));
    allIds = cell(1, numel(worms));
    allIdsFull = cell(1, numel(worms));

    for wi = 1:numel(worms)
        ids = worms(wi).neuron_ids;
        if isempty(ids)
            error('Worm %d is missing neuron_ids; cannot compute shared set.', wi);
        end
        idsFull = normalize_ids_(ids);
        allIdsFull{wi} = idsFull;
        idsIdent = idsFull(is_identified_(idsFull));
        allIds{wi} = idsIdent;
    end

    sharedIds = allIds{1};
    for wi = 2:numel(allIds)
        sharedIds = intersect(sharedIds, allIds{wi}, 'stable');
    end

    for wi = 1:numel(worms)
        [isMember, idx] = ismember(sharedIds, allIdsFull{wi});
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
