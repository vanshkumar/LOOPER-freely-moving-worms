classdef kato_helpers
    % Shared helpers for Kato LOOPER scripts (single + all).
    methods (Static)
        function params = default_kato_params()
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

        function [rawTrain, detrendSpec, split] = prepare_training(worm, trainOnFirstHalf, applyDetrend)
            rawFull = worm.RawData;
            T = size(rawFull, 2);
            [preIdx, postIdx, splitIdx] = looper_helpers.half_split_indices(T);

            if applyDetrend
                if trainOnFirstHalf
                    detrendSpec = looper_helpers.detrend_fit(rawFull, preIdx);
                else
                    detrendSpec = looper_helpers.detrend_fit(rawFull, []);
                end
                rawFull = looper_helpers.detrend_apply(rawFull, detrendSpec);
            else
                detrendSpec = struct('enabled', false, 'mode', 'off');
            end

            if trainOnFirstHalf
                rawTrain = rawFull(:, preIdx);
            else
                rawTrain = rawFull;
            end

            split = struct('enabled', trainOnFirstHalf, ...
                'preIdx', preIdx, 'postIdx', postIdx, 'split_idx', splitIdx);
        end

    end
end
