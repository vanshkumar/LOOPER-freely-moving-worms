classdef looper_eval_helpers
    % Shared helpers for LOOPER eval summaries (Kato + Atanas).
    methods (Static)
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
                'cycles_per_min', nan, 'median_d', nan, 'mean_segment_len', nan);
            if isfield(diag, 'phase_metrics')
                pm = diag.phase_metrics;
                if isfield(pm, 'frac_small'); phase.frac_small = pm.frac_small; end
                if isfield(pm, 'dtheta_var'); phase.dtheta_var = pm.dtheta_var; end
                if isfield(pm, 'cycles_per_min'); phase.cycles_per_min = pm.cycles_per_min; end
                if isfield(pm, 'median_d'); phase.median_d = pm.median_d; end
                if isfield(pm, 'mean_segment_len'); phase.mean_segment_len = pm.mean_segment_len; end
            end

            summary = struct('uid', uid, 'n_neurons', n_neurons, 'T', T, 'dt_sec', dt, ...
                'unique_loops', uniqueLoops, 'loop_switches', loopSwitches, 'segments', segments, ...
                'validation_score_mean', valMean, 'validation_score_std', valStd, ...
                'phase_frac_small', phase.frac_small, 'phase_var', phase.dtheta_var, ...
                'cycles_per_min', phase.cycles_per_min, 'median_d', phase.median_d, ...
                'mean_segment_len', phase.mean_segment_len);
        end
    end
end
