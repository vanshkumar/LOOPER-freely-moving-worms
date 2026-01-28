function saveData = LOOPER(saveData, shouldRun, dynamics, inputs, outputs, parameters)
%LOOPER Run the LOOPER pipeline without the GUI.
%   This mirrors the LOOPER.mlapp "Run all" steps:
%   preprocessing -> diffusion map -> model reduction -> find loops.
%
%   Usage examples:
%     saveData = LOOPER(saveData, true);
%     LOOPER(saveData, true, [], [], [], params);  % update caller's app/saveData
%
%   Inputs:
%     saveData  : struct containing RawData, TrialData, etc.
%     shouldRun : if true, execute the full pipeline
%     dynamics  : optional raw data override (channels x time)
%     inputs    : optional input channels (channels x time)
%     outputs   : optional output channels (channels x time)
%     parameters: struct of overrides (e.g., NearestNeighbors, PreprocessData.*)

    if nargin < 1 || isempty(saveData)
        saveData = struct;
    end
    if nargin < 2 || isempty(shouldRun)
        shouldRun = true;
    end
    if nargin < 3
        dynamics = [];
    end
    if nargin < 4
        inputs = [];
    end
    if nargin < 5
        outputs = [];
    end
    if nargin < 6
        parameters = [];
    end

    % Ensure LOOPER functions are on the path.
    thisDir = fileparts(mfilename('fullpath'));
    looperDir = fullfile(thisDir, 'LOOPER_github_2020', 'Functions');
    if exist(looperDir, 'dir')
        addpath(looperDir);
        addpath(fullfile(looperDir, 'Library'));
        addpath(fullfile(looperDir, 'Library', 'Tubeplot'));
    end

    % Override raw data if provided.
    if ~isempty(dynamics)
        saveData.RawData = dynamics;
    end
    if ~isempty(inputs)
        saveData.Inputs = inputs;
    end
    if ~isempty(outputs)
        saveData.Outputs = outputs;
    end

    % Apply parameter overrides (including PreprocessData sub-struct).
    if exist('parameters', 'var') && ~isempty(parameters)
        fields = fieldnames(parameters);
        for i = 1:length(fields)
            if strcmp(fields{i}, 'PreprocessData')
                if ~isfield(saveData, 'PreprocessData')
                    saveData.PreprocessData = struct;
                end
                pfields = fieldnames(parameters.PreprocessData);
                for j = 1:length(pfields)
                    saveData.PreprocessData.(pfields{j}) = parameters.PreprocessData.(pfields{j});
                end
            else
                saveData.(fields{i}) = parameters.(fields{i});
            end
        end
    end

    % Defaults (match LOOPER.mlapp).
    if ~isfield(saveData, 'PreprocessData') || isempty(saveData.PreprocessData)
        saveData.PreprocessData = struct;
    end
    if ~isfield(saveData.PreprocessData, 'ZScore');      saveData.PreprocessData.ZScore = false; end
    if ~isfield(saveData.PreprocessData, 'Smoothing');   saveData.PreprocessData.Smoothing = 0; end
    if ~isfield(saveData.PreprocessData, 'DelayTime');   saveData.PreprocessData.DelayTime = 0; end
    if ~isfield(saveData.PreprocessData, 'DelayCount');  saveData.PreprocessData.DelayCount = 0; end
    if ~isfield(saveData.PreprocessData, 'InputLambda'); saveData.PreprocessData.InputLambda = 1; end
    if ~isfield(saveData.PreprocessData, 'OutputLambda'); saveData.PreprocessData.OutputLambda = 1; end
    if ~isfield(saveData.PreprocessData, 'TrialLambda'); saveData.PreprocessData.TrialLambda = 1; end
    if ~isfield(saveData.PreprocessData, 'MergeStarts'); saveData.PreprocessData.MergeStarts = false; end
    if ~isfield(saveData.PreprocessData, 'MergeEnds');   saveData.PreprocessData.MergeEnds = false; end

    if ~isfield(saveData, 'NearestNeighbors');   saveData.NearestNeighbors = 6; end
    if ~isfield(saveData, 'UseLocalDimensions'); saveData.UseLocalDimensions = false; end
    if ~isfield(saveData, 'RepopulateDensity');  saveData.RepopulateDensity = 0.95; end
    if ~isfield(saveData, 'MinimumReturnTime');  saveData.MinimumReturnTime = 10; end
    if ~isfield(saveData, 'PutativeClusterCounts'); saveData.PutativeClusterCounts = [100 80 60 50 40 30 20 17 15 12 10 8]; end
    if ~isfield(saveData, 'DistanceType');       saveData.DistanceType = 'correlation'; end
    if ~isfield(saveData, 'MaxCheckTime');       saveData.MaxCheckTime = 10; end
    if ~isfield(saveData, 'PutativeLoopCounts'); saveData.PutativeLoopCounts = [5 4 3 2]; end
    if ~isfield(saveData, 'UseTerminalState');   saveData.UseTerminalState = false; end
    if ~isfield(saveData, 'TotalStates');        saveData.TotalStates = 100; end

    % Provide TrialData if missing (single continuous trial).
    if ~isfield(saveData, 'TrialData') || isempty(saveData.TrialData)
        if isfield(saveData, 'RawData') && ~isempty(saveData.RawData)
            saveData.TrialData = ones(1, size(saveData.RawData, 2));
        else
            saveData.TrialData = [];
        end
    end

    if ~shouldRun
        saveData = assignBack(saveData);
        return;
    end

    % --------------------
    % Step 1: Preprocess
    % --------------------
    rawData = saveData.RawData;
    rawInputs = [];
    rawOutputs = [];
    if isfield(saveData, 'Inputs');  rawInputs = saveData.Inputs;  end
    if isfield(saveData, 'Outputs'); rawOutputs = saveData.Outputs; end

    [tempData, trialData, procesedTrialSwitches, dataMean, dataSTD] = preprocessData( ...
        rawData, rawInputs, saveData.PreprocessData.InputLambda, ...
        rawOutputs, saveData.PreprocessData.OutputLambda, ...
        saveData.TrialData, saveData.PreprocessData.TrialLambda, ...
        true, saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, ...
        saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, [], []);

    saveData.DataMean = dataMean;
    saveData.DataSTD = dataSTD;
    saveData.TrialTimeSeries = trialData;
    saveData.TrialSwitches = procesedTrialSwitches;
    saveData.TimeSeries = tempData;

    % --------------------
    % Step 2: Diffusion map
    % --------------------
    timeSeriesData = saveData.TimeSeries;
    trialTimeSeries = saveData.TrialTimeSeries;
    trialEnds = saveData.TrialSwitches;
    nearestNeighbors = saveData.NearestNeighbors;
    useLocalDimensions = saveData.UseLocalDimensions;
    repopulateDensity = saveData.RepopulateDensity;
    minReturnTime = saveData.MinimumReturnTime;

    cachePath = '';
    if isfield(saveData, 'CacheDiffusionMapPath')
        cachePath = saveData.CacheDiffusionMapPath;
    end

    useCache = false;
    if ~isempty(cachePath) && exist(cachePath, 'file')
        loaded = load(cachePath);
        if isfield(loaded, 'cache')
            loaded = loaded.cache;
        end
        requiredFields = {'asymmetricProbabilities','finalIndicies','finalDynamicsStream', ...
            'stateValidities','diffusionMapIndicies','stateHasNext','pcaBasis'};
        if all(isfield(loaded, requiredFields))
            asymmetricProbabilities = loaded.asymmetricProbabilities;
            finalIndicies = loaded.finalIndicies;
            finalDynamicsStream = loaded.finalDynamicsStream;
            stateValidities = loaded.stateValidities;
            diffusionMapIndicies = loaded.diffusionMapIndicies;
            stateHasNext = loaded.stateHasNext;
            pcaBasis = loaded.pcaBasis;
            useCache = true;
        else
            warning('Diffusion cache missing required fields. Recomputing diffusion map.');
        end
    end

    if ~useCache
        buildDiffusionMap;
        if ~isempty(cachePath)
            cacheDir = fileparts(cachePath);
            if ~isempty(cacheDir) && ~exist(cacheDir, 'dir')
                mkdir(cacheDir);
            end
            cache = struct;
            cache.asymmetricProbabilities = asymmetricProbabilities;
            cache.finalIndicies = finalIndicies;
            cache.finalDynamicsStream = finalDynamicsStream;
            cache.stateValidities = stateValidities;
            cache.diffusionMapIndicies = diffusionMapIndicies;
            cache.stateHasNext = stateHasNext;
            if exist('pcaBasis', 'var')
                cache.pcaBasis = pcaBasis;
            end
            cache.created = datestr(now);
            try
                save(cachePath, 'cache');
            catch saveErr
                warning('Failed to save diffusion cache: %s', saveErr.message);
            end
        end
    end

    saveData.AsymmetricMap = asymmetricProbabilities;
    saveData.FinalIndicies = finalIndicies;
    saveData.FinalStream = finalDynamicsStream;
    saveData.ValidStates = stateValidities;
    saveData.DiffusionMapIndicies = diffusionMapIndicies;
    saveData.StateHasNext = stateHasNext;

    % --------------------
    % Step 3: Model reduction
    % --------------------
    asymmetricProbabilities = saveData.AsymmetricMap;
    finalDynamicsStream = saveData.FinalStream;
    clusterCounts = saveData.PutativeClusterCounts;
    distanceType = saveData.DistanceType;
    maxCheckTime = saveData.MaxCheckTime;
    stateHasNext = saveData.StateHasNext;

    reduceMatrix;

    if ~exist('hadError', 'var') || hadError == 0
        saveData.ReducedMatrix = finalReducedMatrix;
        saveData.ClusterIDs = clusterIDs;
        saveData.ClusterMeansPCA = clusterMeansPCA;
        saveData.ClusterMeans = clusterMeans;
        saveData.CountClusters = countClusters;
    else
        saveData = assignBack(saveData);
        return;
    end

    % --------------------
    % Step 4: Find loops
    % --------------------
    reducedMatrix = saveData.ReducedMatrix;
    clusterIDs = saveData.ClusterIDs;
    stateValidities = saveData.ValidStates;
    clusterMeansPCA = saveData.ClusterMeansPCA;
    clusterMeans = saveData.ClusterMeans;
    countClusters = saveData.CountClusters;
    putativeLoopCounts = saveData.PutativeLoopCounts;
    shouldUseTerminalState = saveData.UseTerminalState;
    totalClusters = saveData.TotalStates;
    trialSwitchTimes = saveData.TrialSwitches;
    selectingLoops = 0;
    app = []; %#ok<NASGU> % placeholder for drawLoops; unused when selectingLoops=0

    if size(saveData.FinalStream,2) > 2
        [pcaBasis, ~] = pca(saveData.FinalStream, 'NumComponents', 3); %#ok<ASGLU>
    else
        pcaBasis = eye(size(saveData.FinalStream,2)); %#ok<NASGU>
    end

    buildMinimalModelFromMatrix;

    saveData.BestStateCount = bestStateCount;
    saveData.LogLikelihoods = allLikelihoods;
    saveData.BestLoopCount = bestLoopCount;
    saveData.BestModel = bestModel;
    saveData.BestEmission = bestEmission;
    saveData.BestLoopAssignments = bestLoopAssignments;
    saveData.BestStateMap = bestStateMap;

    saveData = assignBack(saveData);
end

function saveData = assignBack(saveData)
%ASSIGNBACK Update caller workspace for compatibility with legacy scripts.
    try
        if evalin('caller', 'exist(''app'', ''var'')')
            app = evalin('caller', 'app');
            if isstruct(app) || isobject(app)
                app.SavedData = saveData;
            else
                app = struct('SavedData', saveData);
            end
            assignin('caller', 'app', app);
        else
            assignin('caller', 'app', struct('SavedData', saveData));
        end
        assignin('caller', 'saveData', saveData);
    catch
        % Best effort only; ignore if caller workspace is not writable.
    end
end
