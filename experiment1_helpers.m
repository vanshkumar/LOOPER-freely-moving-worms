classdef experiment1_helpers
    % Shared helpers for Experiment 1 scripts.
    methods (Static)
        function [preIdx, postIdx] = normalize_ranges(rangesRaw)
            % Accepts cell or numeric ranges and returns row vectors.
            if iscell(rangesRaw)
                preIdx = rangesRaw{1};
                postIdx = rangesRaw{2};
                return;
            end

            % If numeric, assume rows or columns.
            if size(rangesRaw, 1) == 2
                preIdx = rangesRaw(1, :);
                postIdx = rangesRaw(2, :);
            else
                preIdx = rangesRaw(:, 1);
                postIdx = rangesRaw(:, 2);
            end
            preIdx = preIdx(:)';
            postIdx = postIdx(:)';
        end

        function [alpha, theta, stateIdx, d] = project_to_model(rawData, saveData)
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

        function dtheta = circ_diff_bins(thetaSeq, nBins)
            if numel(thetaSeq) < 2
                dtheta = nan(size(thetaSeq));
                return;
            end
            deltas = diff(thetaSeq);
            dtheta = experiment1_helpers.wrap_diff_bins(deltas, nBins);
        end

        function d = wrap_diff_bins(d, nBins)
            % Wrap to [-nBins/2, nBins/2).
            d = mod(d + nBins/2, nBins) - nBins/2;
        end

        function setup_figs()
            % Use docked tabs for consistent LOOPER plotting behavior.
            set(0, 'DefaultFigureWindowStyle', 'docked');
            set(0, 'DefaultFigureVisible', 'on');
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

        function params = default_looper_params_atanas()
            % Paper-style C. elegans settings, scaled to Atanas dt (~0.6s).
            params = [];
            params.NearestNeighbors = 8;
            params.UseLocalDimensions = true;
            params.RepopulateDensity = 0.95;
            params.MinimumReturnTime = 4;
            params.DistanceType = 'correlation';
            params.MaxCheckTime = 4;
            params.TotalStates = 25;
            params.UseTerminalState = false;
            % Wider loop-count range to avoid edge minima (per paper note).
            params.PutativeLoopCounts = [8 7 6 5 4 3 2];
            params.PreprocessData.ZScore = false;
            params.PreprocessData.Smoothing = 1;
            params.PreprocessData.DelayTime = 4;
            params.PreprocessData.DelayCount = 5;
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

        function warn_ranges_contiguous(preIdx, postIdx, uid)
            if nargin < 3
                uid = '';
            end
            if numel(preIdx) <= 50 || numel(postIdx) <= 50
                warning('ranges_raw looks like start/end, not index vectors (%s).', uid);
            end
            if preIdx(end) + 1 ~= postIdx(1)
                warning('ranges_raw split not contiguous for %s.', uid);
            end
        end

        function metrics = phase_continuity_metrics(alpha, theta, d, nBins, dt)
            % Summarize phase continuity within constant-alpha segments.
            if nargin < 5 || isempty(dt)
                dt = nan;
            end
            if isempty(alpha) || isempty(theta)
                metrics = struct('frac_small', nan, 'dtheta_var', nan, ...
                    'cycles_per_min', nan, 'median_d', nan, ...
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
                dtheta = experiment1_helpers.wrap_diff_bins(diff(thetaSeg), nBins);
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

        function summary = summary_metrics(alpha, theta, d, relTime, delta_d, worm, t_on, dt, nBins)
            % Build summary metrics for Experiment 1 evals.
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

            phaseMetrics = experiment1_helpers.phase_continuity_metrics(alpha, theta, d, nBins, dt);
            diag = struct('phase_metrics', phaseMetrics);
            meta = struct('uid', worm.uid, 'n_neurons', worm.num_neurons, 'T', worm.T, 'dt_sec', dt);

            uniqueLoops = numel(unique(alpha));
            loopSwitches = sum(diff(alpha) ~= 0);
            segments = 1 + loopSwitches;

            summary = looper_eval_helpers.base_summary([], diag, meta);
            summary.unique_loops = uniqueLoops;
            summary.loop_switches = loopSwitches;
            summary.segments = segments;
            summary.split_idx = t_on;
            summary.mean_pre = meanPre;
            summary.mean_post = meanPost;
            summary.post_slope = postSlope;
            summary.d_peak = dPeak;
        end
    end
end
