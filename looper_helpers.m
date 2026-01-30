classdef looper_helpers
    % Shared helpers for LOOPER eval + projection utilities.
    methods (Static)
        function [preIdx, postIdx, splitIdx] = half_split_indices(T)
            % Half-split indices for stationarity runs (shared across datasets).
            splitIdx = floor(T / 2);
            preIdx = 1:splitIdx;
            postIdx = (splitIdx + 1):T;
        end

        function summary = base_summary(saveData, diag, meta)
            % Build summary fields shared across evals.
            if nargin < 3
                meta = struct;
            end

            uid = '';
            if isfield(meta, 'uid') && ~isempty(meta.uid)
                uid = meta.uid;
            end

            n_neurons = nan;
            T = nan;
            if isfield(meta, 'n_neurons') && ~isempty(meta.n_neurons)
                n_neurons = meta.n_neurons;
            elseif isfield(saveData, 'RawData') && ~isempty(saveData.RawData)
                n_neurons = size(saveData.RawData, 1);
            end
            if isfield(meta, 'T') && ~isempty(meta.T)
                T = meta.T;
            elseif isfield(saveData, 'RawData') && ~isempty(saveData.RawData)
                T = size(saveData.RawData, 2);
            end

            dt = nan;
            if isfield(meta, 'dt_sec') && ~isempty(meta.dt_sec)
                dt = meta.dt_sec;
            end

            uniqueLoops = nan;
            loopSwitches = nan;
            segments = nan;
            if isfield(saveData, 'BestStateMap') && ~isempty(saveData.BestStateMap)
                alpha = saveData.BestStateMap(:,1);
                uniqueLoops = numel(unique(alpha));
                loopSwitches = sum(diff(alpha) ~= 0);
                segments = 1 + loopSwitches;
            end

            valMean = nan;
            valStd = nan;
            if isfield(diag, 'validationScoreMean')
                valMean = diag.validationScoreMean;
            end
            if isfield(diag, 'validationScoreStd')
                valStd = diag.validationScoreStd;
            end

            phase = struct('frac_small', nan, 'dtheta_var', nan, ...
                'phase_speed_bins_per_min', nan, 'median_d', nan, 'mean_segment_len', nan);
            if isfield(diag, 'phase_metrics')
                pm = diag.phase_metrics;
                if isfield(pm, 'frac_small'); phase.frac_small = pm.frac_small; end
                if isfield(pm, 'dtheta_var'); phase.dtheta_var = pm.dtheta_var; end
                if isfield(pm, 'phase_speed_bins_per_min'); phase.phase_speed_bins_per_min = pm.phase_speed_bins_per_min; end
                if isfield(pm, 'median_d'); phase.median_d = pm.median_d; end
                if isfield(pm, 'mean_segment_len'); phase.mean_segment_len = pm.mean_segment_len; end
            end

            reconCorr = nan;
            if isfield(diag, 'recon_corr_full')
                reconCorr = diag.recon_corr_full;
            end

            summary = struct('uid', uid, 'n_neurons', n_neurons, 'T', T, 'dt_sec', dt, ...
                'unique_loops', uniqueLoops, 'loop_switches', loopSwitches, 'segments', segments, ...
                'validation_score_mean', valMean, 'validation_score_std', valStd, ...
                'recon_corr_full', reconCorr, ...
                'phase_frac_small', phase.frac_small, 'phase_var', phase.dtheta_var, ...
                'phase_speed_bins_per_min', phase.phase_speed_bins_per_min, 'median_d', phase.median_d, ...
                'mean_segment_len', phase.mean_segment_len, ...
                'split_idx', nan, 'mean_pre', nan, 'mean_post', nan, ...
                'post_slope', nan, 'd_peak', nan, 'recon_corr_post', nan);
        end

        function [alpha, theta, stateIdx, d, X, clusterMeans] = project_to_model(rawData, saveData)
            % Project full data onto a learned LOOPER scaffold via emission matching.
            rawTrialData = ones(1, size(rawData, 2));
            [tempData, ~, ~] = preprocessData( ...
                rawData, [], 0, [], 0, rawTrialData, 0, false, ...
                saveData.PreprocessData.Smoothing, saveData.PreprocessData.ZScore, ...
                saveData.PreprocessData.DelayTime, saveData.PreprocessData.DelayCount, ...
                saveData.DataMean, saveData.DataSTD);

            X = tempData';

            finalStream = saveData.FinalStream;
            if size(finalStream,1) ~= size(saveData.BestStateMap,1)
                warning('FinalStream and BestStateMap lengths differ; projection may be inconsistent.');
            end
            numStates = size(saveData.BestLoopAssignments, 1);
            clusterMeans = zeros(numStates, size(finalStream, 2));
            clusterStds = zeros(numStates, size(finalStream, 2));
            globalMean = mean(finalStream, 1);
            globalStd = std(finalStream, [], 1);
            globalStd(~isfinite(globalStd) | globalStd == 0) = 1;

            for i = 1:numStates
                thisLoopPosition = saveData.BestLoopAssignments(i,:);
                thisIdx = find(ismember(saveData.BestStateMap, thisLoopPosition, 'rows'));
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

            if size(clusterMeans,2) ~= size(X,2)
                error('Projection mismatch: BestEmission has %d dims, but processed data has %d.', ...
                    size(clusterMeans,2), size(X,2));
            end

            trim = saveData.PreprocessData.DelayCount * saveData.PreprocessData.DelayTime;
            expectedLen = size(rawData,2) - trim;
            if abs(size(X,1) - expectedLen) > 1
                warning('Projection length (%d) differs from expected (%d) by >1. Check delay embedding.', ...
                    size(X,1), expectedLen);
            end

            dists = zeros(size(clusterMeans,1), size(X,1));
            for i = 1:size(clusterMeans,1)
                denom = clusterStds(i,:);
                denom(denom == 0) = 1;
                z = (X - clusterMeans(i,:)) ./ denom;
                dists(i,:) = sum(z.^2, 2)';
            end
            [~, stateIdx] = min(dists, [], 1);

            alpha = saveData.BestLoopAssignments(stateIdx, 1);
            theta = saveData.BestLoopAssignments(stateIdx, 2);

            % Use the same z-scored distance used for assignment.
            d = sqrt(min(dists, [], 1))';
        end

        function d = wrap_diff_bins(d, nBins)
            % Wrap to [-nBins/2, nBins/2).
            d = mod(d + nBins/2, nBins) - nBins/2;
        end

        function detrend = detrend_fit(rawData, fitIdx)
            % Fit per-neuron linear trend (fit on fitIdx, apply to all).
            if nargin < 2 || isempty(fitIdx)
                fitIdx = 1:size(rawData, 2);
            end
            t = 1:size(rawData, 2);
            tFit = t(fitIdx);
            slopes = zeros(size(rawData, 1), 1);
            intercepts = zeros(size(rawData, 1), 1);
            for i = 1:size(rawData, 1)
                y = rawData(i, fitIdx);
                p = polyfit(tFit, y, 1);
                slopes(i) = p(1);
                intercepts(i) = p(2);
            end
            detrend = struct('enabled', true, 'mode', 'fit_pre', ...
                'slopes', slopes, 'intercepts', intercepts);
        end

        function rawDet = detrend_apply(rawData, detrend)
            % Apply saved linear detrend to rawData.
            t = 1:size(rawData, 2);
            trend = detrend.slopes * t + detrend.intercepts;
            rawDet = rawData - trend;
        end

        function out = compute_delta_d(d, dt, saveData, t_on, preIdx, postIdx)
            % Compute baseline-normalized distance around an event window.
            trim = saveData.PreprocessData.DelayCount * saveData.PreprocessData.DelayTime;
            procLen = numel(d);
            tOnProc = t_on - trim;
            if tOnProc < 1 || tOnProc > procLen
                error(['Event index falls outside delay-embedded window. ', ...
                    'Reduce DelayCount/DelayTime or adjust indices.']);
            end

            preIdxProc = preIdx(preIdx > trim) - trim;
            postIdxProc = postIdx(postIdx > trim) - trim;
            preIdxProc = preIdxProc(preIdxProc >= 1 & preIdxProc <= procLen);
            postIdxProc = postIdxProc(postIdxProc >= 1 & postIdxProc <= procLen);
            preIdxProc = preIdxProc(preIdxProc < tOnProc);
            if isempty(preIdxProc)
                preIdxProc = 1:max(1, tOnProc - 1);
            end
            if isempty(postIdxProc)
                postIdxProc = tOnProc:procLen;
            end

            idxStart = preIdxProc(1);
            idxEnd = postIdxProc(end);
            relIdx = idxStart:idxEnd;
            relTime = (relIdx - tOnProc) * dt;
            baseline = median(d(preIdxProc), 'omitnan');
            delta_d = d(relIdx) - baseline;

            out = struct;
            out.trim = trim;
            out.tOnProc = tOnProc;
            out.preIdxProc = preIdxProc;
            out.postIdxProc = postIdxProc;
            out.relIdx = relIdx;
            out.relTime = relTime;
            out.baseline = baseline;
            out.delta_d = delta_d;
        end

        function metrics = phase_continuity_metrics(alpha, theta, d, nBins, dt)
            % Summarize phase continuity within constant-alpha segments.
            if nargin < 5 || isempty(dt)
                dt = nan;
            end
            if isempty(alpha) || isempty(theta)
                metrics = struct('frac_small', nan, 'dtheta_var', nan, ...
                    'phase_speed_bins_per_min', nan, 'median_d', nan, ...
                    'segments', 0, 'mean_segment_len', nan, 'nBins', nBins);
                return;
            end

            alpha = alpha(:);
            theta = theta(:);
            if ~isempty(d)
                d = d(:);
            else
                d = nan(size(alpha));
            end

            changeIdx = find(diff(alpha) ~= 0);
            edges = [1; changeIdx + 1; numel(alpha) + 1];

            dthetaAll = [];
            dAll = [];
            totalCycles = 0;
            totalPhaseSteps = 0;
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
                dtheta = looper_helpers.wrap_diff_bins(diff(thetaSeg), nBins);
                dthetaAll = [dthetaAll; dtheta(:)]; %#ok<AGROW>
                totalPhaseSteps = totalPhaseSteps + sum(abs(dtheta));

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
                metrics.phase_speed_bins_per_min = nan;
            else
                metrics.phase_speed_bins_per_min = totalPhaseSteps / (totalTimeSec / 60);
            end
            metrics.median_d = median(dAll, 'omitnan');
            if isempty(segLens)
                metrics.mean_segment_len = nan;
            else
                metrics.mean_segment_len = mean(segLens);
            end
        end

        function summary = summary_metrics(alpha, theta, d, relTime, delta_d, wormMeta, ...
                t_on, dt, nBins, diagIn, reconCorrPost)
            % Build summary metrics for split-half evals.
            if nargin < 10 || isempty(diagIn)
                diagIn = struct;
            end
            if nargin < 11 || isempty(reconCorrPost)
                reconCorrPost = nan;
            end
            preMask = relTime < 0;
            postMask = relTime >= 0;
            meanPre = mean(delta_d(preMask), 'omitnan');
            meanPost = mean(delta_d(postMask), 'omitnan');
            dPeak = max(delta_d(postMask), [], 'omitnan');

            postSlope = nan;
            if any(postMask)
                p = polyfit(relTime(postMask), delta_d(postMask), 1);
                postSlope = p(1);
            end

            phaseMetrics = looper_helpers.phase_continuity_metrics(alpha, theta, d, nBins, dt);
            diag = struct('phase_metrics', phaseMetrics);
            if isfield(diagIn, 'validationScoreMean'); diag.validationScoreMean = diagIn.validationScoreMean; end
            if isfield(diagIn, 'validationScoreStd'); diag.validationScoreStd = diagIn.validationScoreStd; end
            if isfield(diagIn, 'recon_corr_full'); diag.recon_corr_full = diagIn.recon_corr_full; end

            meta = struct('uid', '', 'n_neurons', nan, 'T', nan, 'dt_sec', dt);
            if isstruct(wormMeta)
                if isfield(wormMeta, 'uid'); meta.uid = wormMeta.uid; end
                if isfield(wormMeta, 'num_neurons'); meta.n_neurons = wormMeta.num_neurons; end
                if isfield(wormMeta, 'n_neurons'); meta.n_neurons = wormMeta.n_neurons; end
                if isfield(wormMeta, 'T'); meta.T = wormMeta.T; end
            end

            uniqueLoops = numel(unique(alpha));
            loopSwitches = sum(diff(alpha) ~= 0);
            segments = 1 + loopSwitches;

            summary = looper_helpers.base_summary([], diag, meta);
            summary.unique_loops = uniqueLoops;
            summary.loop_switches = loopSwitches;
            summary.segments = segments;
            summary.split_idx = t_on;
            summary.mean_pre = meanPre;
            summary.mean_post = meanPost;
            summary.post_slope = postSlope;
            summary.d_peak = dPeak;
            summary.recon_corr_post = reconCorrPost;
        end

        function evalOut = eval_stationarity(rawFull, saveData, preIdx, postIdx, dt, wormMeta, diag)
            % Shared stationarity evaluation helper.
            if nargin < 6 || isempty(wormMeta)
                wormMeta = struct;
            end
            if nargin < 7
                diag = struct;
            end
            t_on = postIdx(1);

            if isfield(saveData, 'Detrend') && isfield(saveData.Detrend, 'enabled') && saveData.Detrend.enabled
                rawFull = looper_helpers.detrend_apply(rawFull, saveData.Detrend);
            end

            [alpha, theta, stateIdx, d, X, clusterMeans] = looper_helpers.project_to_model(rawFull, saveData);

            if isempty(dt) || ~isfinite(dt)
                dt = 1;
            end

            out = looper_helpers.compute_delta_d(d, dt, saveData, t_on, preIdx, postIdx);

            nBins = max(saveData.BestLoopAssignments(:,2));
            if isempty(nBins) || ~isfinite(nBins)
                nBins = max(theta);
            end

            reconCorrPost = nan;
            if ~isempty(clusterMeans) && ~isempty(stateIdx) && ~isempty(X)
                procStream = X;
                reconStream = clusterMeans(stateIdx, :);
                postIdxProc = out.postIdxProc;
                postIdxProc = postIdxProc(postIdxProc >= 1 & postIdxProc <= size(procStream, 1));
                if ~isempty(postIdxProc) && size(procStream, 2) == size(reconStream, 2)
                    postFinal = procStream(postIdxProc, :);
                    postRecon = reconStream(postIdxProc, :);
                    reconCorrPost = corr(postFinal(:), postRecon(:));
                end
            end

            summary = looper_helpers.summary_metrics( ...
                alpha, theta, d, out.relTime, out.delta_d, wormMeta, ...
                t_on, dt, nBins, diag, reconCorrPost);

            evalOut = struct;
            evalOut.alpha = alpha;
            evalOut.theta = theta;
            evalOut.stateIdx = stateIdx;
            evalOut.d = d;
            evalOut.relTime = out.relTime;
            evalOut.delta_d = out.delta_d;
            evalOut.t_on = t_on;
            evalOut.tOnProc = out.tOnProc;
            evalOut.preIdx = preIdx;
            evalOut.postIdx = postIdx;
            evalOut.dt = dt;
            evalOut.nBins = nBins;
            evalOut.summary = summary;
            evalOut.postIdxProc = out.postIdxProc;
            evalOut.preIdxProc = out.preIdxProc;
        end

        function save_stationarity_artifacts(evalOut, worm, plotDir, projDir, opts)
            % Save recovery/loop plots and projection mat for stationarity runs.
            if nargin < 5
                opts = struct;
            end
            if ~isfield(opts, 'tag'); opts.tag = ''; end
            if ~isfield(opts, 'makeLoopPlot'); opts.makeLoopPlot = true; end

            if isempty(plotDir)
                plotDir = pwd;
            end
            if nargin < 4 || isempty(projDir)
                projDir = plotDir;
            end
            if ~exist(plotDir, 'dir')
                mkdir(plotDir);
            end
            if ~exist(projDir, 'dir')
                mkdir(projDir);
            end

            uid = '';
            if isstruct(worm) && isfield(worm, 'uid')
                uid = worm.uid;
            end

            tag = opts.tag;
            if isempty(tag)
                tag = uid;
            end
            safeTag = regexprep(tag, '[^A-Za-z0-9_-]', '_');

            recoveryName = sprintf('%s_stationarity_recovery.png', safeTag);
            loopName = sprintf('%s_loop_assignment.png', safeTag);
            projName = sprintf('%s_stationarity_projection.mat', safeTag);

            figRecovery = figure('Name', sprintf('Stationarity recovery (%s)', uid), ...
                'NumberTitle', 'off');
            plot(evalOut.relTime, evalOut.delta_d, 'LineWidth', 1.5);
            xline(0, '--k');
            xlabel('Seconds from split');
            ylabel('\\Delta d (baseline-normalized)');
            title(sprintf('Stationarity recovery (worm %s)', uid));
            saveas(figRecovery, fullfile(plotDir, recoveryName));

            if opts.makeLoopPlot
                figLoops = figure('Name', sprintf('Loop assignment (%s)', uid), ...
                    'NumberTitle', 'off');
                plot(evalOut.alpha, 'LineWidth', 1.2);
                xline(evalOut.tOnProc, '--k');
                xlabel('Time (embedded frames)');
                ylabel('Loop assignment');
                title(sprintf('Loop assignment over time (worm %s)', uid));
                saveas(figLoops, fullfile(plotDir, loopName));
            end

            alpha = evalOut.alpha; %#ok<NASGU>
            theta = evalOut.theta; %#ok<NASGU>
            stateIdx = evalOut.stateIdx; %#ok<NASGU>
            d = evalOut.d; %#ok<NASGU>
            save(fullfile(projDir, projName), 'alpha', 'theta', 'stateIdx', 'd', 'evalOut', 'worm');
        end
    end
end
