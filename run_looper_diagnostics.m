function diag = run_looper_diagnostics(saveData, opts)
%RUN_LOOPER_DIAGNOSTICS Run LOOPER's built-in diagnostics on saveData.
%
% Usage:
%   diag = run_looper_diagnostics(saveData);
%   diag = run_looper_diagnostics('results/experiment1_single/xxxx.mat');
%   diag = run_looper_diagnostics(saveData, struct('preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));
%
% Options (opts fields, all optional):
%   - outDir: output folder for figures + diag mat
%   - tag: label for filenames
%   - preIdxRaw/postIdxRaw/tOnRaw: raw indices (pre-embedding)
%   - preIdxProc/postIdxProc/tOnProc: indices in embedded timebase
%   - saveFigs: true/false (default true)
%   - checkTrainMatch: true/false (default true)

    if nargin < 2
        opts = struct;
    end

    if ischar(saveData) || isstring(saveData)
        loaded = load(saveData);
        if isfield(loaded, 'saveData')
            saveData = loaded.saveData;
        else
            error('MAT file must contain variable saveData.');
        end
    end

    addpath(fullfile(pwd, 'LOOPER_github_2020', 'Functions'));
    addpath(fullfile(pwd, 'LOOPER_github_2020', 'Functions', 'Library'));
    addpath(fullfile(pwd, 'LOOPER_github_2020', 'Functions', 'Library', 'Tubeplot'));
    addpath(fullfile(pwd, 'LOOPER_OSF_hosted_2022'));

    if ~isfield(opts, 'saveFigs')
        opts.saveFigs = true;
    end
    if ~isfield(opts, 'cloneFigs')
        opts.cloneFigs = true;
    end
    if ~isfield(opts, 'closeOriginalFigs')
        opts.closeOriginalFigs = true;
    end
    if ~isfield(opts, 'tag') || isempty(opts.tag)
        opts.tag = 'looper';
    end
    if ~isfield(opts, 'outDir') || isempty(opts.outDir)
        opts.outDir = fullfile(pwd, 'results', 'diagnostics', opts.tag);
    end
    if ~exist(opts.outDir, 'dir')
        mkdir(opts.outDir);
    end
    if ~isfield(opts, 'checkTrainMatch')
        opts.checkTrainMatch = true;
    end

    numTrial = 1;
    if isfield(saveData, 'TrialData') && ~isempty(saveData.TrialData)
        numTrial = max(saveData.TrialData);
    end

    if ~isfield(saveData, 'RawTrialSwitches') || isempty(saveData.RawTrialSwitches)
        if numTrial == 1
            saveData.RawTrialSwitches = [];
        elseif isfield(saveData, 'TrialSwitches') && ~isempty(saveData.TrialSwitches)
            reducedSize = (size(saveData.TrialData,2) - size(saveData.FinalStream,1)) / numTrial;
            saveData.RawTrialSwitches = saveData.TrialSwitches + (0:numTrial-1) * reducedSize;
        else
            saveData.RawTrialSwitches = [];
        end
    end

    if ~isempty(saveData.RawTrialSwitches)
        maxLen = size(saveData.BestStateMap, 1);
        saveData.RawTrialSwitches = saveData.RawTrialSwitches(saveData.RawTrialSwitches < maxLen);
    end

    diag = struct;
    diag.tag = opts.tag;

    % --- Plot loops (tube visualization) ---
    plotLoops;
    if opts.cloneFigs
        figHandle = findobj('Type', 'figure', 'Number', 10000);
        clone_fig_(figHandle, sprintf('%s_plotLoops', opts.tag), opts.closeOriginalFigs);
    end
    if opts.saveFigs
        saveas(gcf, fullfile(opts.outDir, sprintf('%s_plotLoops.png', opts.tag)));
    end

    % --- Asymmetric transition matrix ---
    if isfield(saveData, 'AsymmetricMap') && ~isempty(saveData.AsymmetricMap)
        figTrans = figure(10006); clf;
        imagesc(saveData.AsymmetricMap);
        axis square;
        title('Asymmetric transition matrix');
        xlabel('To state');
        ylabel('From state');
        colormap(parula(256));
        if opts.cloneFigs
            clone_fig_(figTrans, sprintf('%s_transition_matrix', opts.tag), opts.closeOriginalFigs);
        end
        if opts.saveFigs
            saveas(gcf, fullfile(opts.outDir, sprintf('%s_transition_matrix.png', opts.tag)));
        end
    end

    % --- Reconstruction R^2 ---
    clear reconstructionPlotIndices;
    plotReconstruction;
    if opts.cloneFigs
        figHandle = findobj('Type', 'figure', 'Number', 1001);
        clone_fig_(figHandle, sprintf('%s_reconstruction', opts.tag), opts.closeOriginalFigs);
    end
    if exist('Rsquared', 'var')
        diag.R2_full = Rsquared;
    end
    if opts.saveFigs
        saveas(gcf, fullfile(opts.outDir, sprintf('%s_reconstruction_full.png', opts.tag)));
    end

    clear reconstructionPlotIndices;

    % --- Validate model (LOOPER internal) ---
    validationModel = saveData.BestModel;
    validationEmission = saveData.BestEmission;
    finalDynamicsStream = saveData.FinalStream;
    originalDataSize = size(saveData.RawData);
    validateModel;
    if opts.cloneFigs
        figHandle = findobj('Type', 'figure', 'Number', 301);
        clone_fig_(figHandle, sprintf('%s_validateModel', opts.tag), opts.closeOriginalFigs);
    end
    if exist('scoreMean', 'var')
        diag.validationScoreMean = scoreMean;
        diag.validationScoreStd = scoreSTD;
    end
    if opts.saveFigs
        saveas(gcf, fullfile(opts.outDir, sprintf('%s_validateModel.png', opts.tag)));
    end

    % --- Loop/phase vs time (OSF script) ---
    if size(saveData.BestStateMap,1) == size(saveData.FinalStream,1)
        trialLengths = [];
        if isfield(saveData, 'RawTrialSwitches') && ~isempty(saveData.RawTrialSwitches)
            trialLengths = diff([0, saveData.RawTrialSwitches, size(saveData.BestStateMap,1)]);
        end

        if ~isempty(trialLengths) && numel(unique(trialLengths)) > 1
            diag.validateDataSkipped = true;
            diag.validateDataReason = 'Variable trial lengths; OSF validateData assumes equal-length trials.';
        else
            SHOULD_VALIDATE = 0;
            totalTrials = numTrial;
            validateData;
            if opts.cloneFigs
                figHandle = findobj('Type', 'figure', 'Number', 1);
                clone_fig_(figHandle, sprintf('%s_loop_phase_time', opts.tag), opts.closeOriginalFigs);
            end
            if opts.saveFigs
                saveas(gcf, fullfile(opts.outDir, sprintf('%s_loop_phase_time.png', opts.tag)));
            end
        end
    else
        diag.validateDataSkipped = true;
        diag.validateDataReason = 'BestStateMap length does not match FinalStream length.';
    end

    % --- Training round-trip check (projection vs BestStateMap) ---
    if opts.checkTrainMatch && isfield(saveData, 'RawData') && ~isempty(saveData.RawData) ...
            && isfield(opts, 'preIdxRaw') && ~isempty(opts.preIdxRaw)
        try
            rawTrain = saveData.RawData(:, opts.preIdxRaw);
            [alphaTrain, thetaTrain] = experiment1_helpers.project_to_model(rawTrain, saveData);
            alphaLooper = saveData.BestStateMap(:,1);
            thetaLooper = saveData.BestStateMap(:,2);
            L = min(numel(alphaTrain), numel(alphaLooper));
            if L > 0
                diag.trainMatchAlpha = mean(alphaTrain(1:L) == alphaLooper(1:L));
                diag.trainMatchTheta = mean(thetaTrain(1:L) == thetaLooper(1:L));
                fprintf('%s train round-trip: alpha=%.3f theta=%.3f\n', ...
                    opts.tag, diag.trainMatchAlpha, diag.trainMatchTheta);
            else
                diag.trainMatchAlpha = nan;
                diag.trainMatchTheta = nan;
            end
        catch err
            diag.trainMatchError = err.message;
            warning('Train round-trip failed for %s: %s', opts.tag, err.message);
        end
    end

    % --- Phase continuity diagnostics (loop-ishness) ---
    if size(saveData.BestStateMap,1) == size(saveData.FinalStream,1)
        alpha = saveData.BestStateMap(:,1);
        theta = saveData.BestStateMap(:,2);
        nBins = max(saveData.BestLoopAssignments(:,2));
        if isempty(nBins) || ~isfinite(nBins)
            nBins = max(theta);
        end

        % Compute per-state mean/std from the training stream (OSF-style).
        stateIdx = ismember_rows_(saveData.BestStateMap, saveData.BestLoopAssignments);
        finalStream = saveData.FinalStream;
        numStates = size(saveData.BestLoopAssignments, 1);
        clusterMeans = zeros(numStates, size(finalStream, 2));
        clusterStds = zeros(numStates, size(finalStream, 2));
        globalMean = mean(finalStream, 1);
        globalStd = std(finalStream, [], 1);
        globalStd(~isfinite(globalStd) | globalStd == 0) = 1;
        for i = 1:numStates
            thisIdx = find(stateIdx == i);
            if isempty(thisIdx)
                clusterMeans(i,:) = globalMean;
                clusterStds(i,:) = globalStd;
            else
                clusterMeans(i,:) = mean(finalStream(thisIdx,:), 1);
                thisStd = std(finalStream(thisIdx,:), [], 1);
                thisStd(~isfinite(thisStd) | thisStd == 0) = globalStd(~isfinite(thisStd) | thisStd == 0);
                clusterStds(i,:) = thisStd;
            end
        end

        d = nan(size(alpha));
        for i = 1:numel(alpha)
            s = stateIdx(i);
            if s < 1
                continue;
            end
            z = (finalStream(i,:) - clusterMeans(s,:)) ./ clusterStds(s,:);
            d(i) = sqrt(sum(z.^2));
        end

        dt = nan;
        if isfield(opts, 'dt_sec') && ~isempty(opts.dt_sec)
            dt = opts.dt_sec;
        end
        diag.phase_metrics = phase_metrics_(alpha, theta, d, nBins, dt);
    else
        diag.phase_metrics = struct('frac_small', nan, 'dtheta_var', nan, ...
            'cycles_per_min', nan, 'median_d', nan, ...
            'segments', 0, 'mean_segment_len', nan, 'nBins', nan);
    end

    % Save diagnostics struct.
    save(fullfile(opts.outDir, sprintf('%s_diag.mat', opts.tag)), 'diag');
end

function clone_fig_(figHandle, figName, closeOriginal)
    if isempty(figHandle)
        return;
    end
    if numel(figHandle) > 1
        figHandle = figHandle(1);
    end
    if ~ishandle(figHandle)
        return;
    end
    newFig = figure('Name', figName, 'NumberTitle', 'off');
    copyobj(get(figHandle, 'Children'), newFig);
    try
        colormap(newFig, colormap(figHandle));
    catch
    end
    if closeOriginal
        close(figHandle);
    end
end

function idx = ismember_rows_(A, B)
    % Returns row indices of A in B (0 if not found).
    [~, idx] = ismember(A, B, 'rows');
end

function metrics = phase_metrics_(alpha, theta, d, nBins, dt)
    if nargin < 5 || isempty(dt)
        dt = nan;
    end
    alpha = alpha(:);
    theta = theta(:);
    d = d(:);

    changeIdx = find(diff(alpha) ~= 0);
    edges = [1; changeIdx + 1; numel(alpha) + 1];

    dthetaAll = [];
    dAll = [];
    totalCycles = 0;
    totalTimeSec = 0;
    segLens = [];

    for i = 1:numel(edges)-1
        segStart = edges(i);
        segEnd = edges(i+1) - 1;
        segLen = segEnd - segStart + 1;
        if segLen < 2
            continue;
        end
        segLens(end+1,1) = segLen; %#ok<AGROW>
        thetaSeg = theta(segStart:segEnd);
        dtheta = wrap_diff_bins_(diff(thetaSeg), nBins);
        dthetaAll = [dthetaAll; dtheta(:)]; %#ok<AGROW>

        if all(isfinite(thetaSeg))
            thetaUnwrap = cumsum([thetaSeg(1); dtheta(:)]);
            cycles = (thetaUnwrap(end) - thetaUnwrap(1)) / nBins;
            totalCycles = totalCycles + abs(cycles);
        end
        totalTimeSec = totalTimeSec + segLen * dt;

        dAll = [dAll; d(segStart:segEnd)]; %#ok<AGROW>
    end

    metrics = struct;
    metrics.nBins = nBins;
    metrics.segments = numel(segLens);
    if isempty(dthetaAll)
        metrics.frac_small = nan;
        metrics.dtheta_var = nan;
    else
        metrics.frac_small = mean(abs(dthetaAll) <= 1, 'omitnan');
        metrics.dtheta_var = var(dthetaAll, 'omitnan');
    end
    if ~isfinite(totalTimeSec) || totalTimeSec <= 0
        metrics.cycles_per_min = nan;
    else
        metrics.cycles_per_min = totalCycles / (totalTimeSec / 60);
    end
    metrics.median_d = median(dAll, 'omitnan');
    if isempty(segLens)
        metrics.mean_segment_len = nan;
    else
        metrics.mean_segment_len = mean(segLens);
    end
end

function d = wrap_diff_bins_(d, nBins)
    d = mod(d + nBins/2, nBins) - nBins/2;
end
