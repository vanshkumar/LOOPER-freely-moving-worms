classdef atanas_helpers
    % Shared helpers for Atanas scripts.
    methods (Static)
        function [rawTrain, detrendSpec, split] = prepare_training(worm, fitFirstHalf, applyDetrend)
            % Prepare LOOPER training data for stationarity split runs.
            rawFull = worm.RawData;
            split = struct('enabled', fitFirstHalf, 'preIdx', [], ...
                'postIdx', [], 'split_idx', []);

            preIdx = [];
            postIdx = [];
            if fitFirstHalf
                [preIdx, postIdx, splitIdx] = looper_helpers.half_split_indices(size(rawFull, 2));
                split.preIdx = preIdx;
                split.postIdx = postIdx;
                split.split_idx = splitIdx;
            end

            if applyDetrend
                if fitFirstHalf
                    detrendSpec = looper_helpers.detrend_fit(rawFull, preIdx);
                else
                    detrendSpec = looper_helpers.detrend_fit(rawFull, []);
                end
                rawFull = looper_helpers.detrend_apply(rawFull, detrendSpec);
            else
                detrendSpec = struct('enabled', false, 'mode', 'off');
            end

            if fitFirstHalf
                rawTrain = rawFull(:, preIdx);
            else
                rawTrain = rawFull;
            end
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

    end
end
