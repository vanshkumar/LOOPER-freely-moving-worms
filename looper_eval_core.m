function looper_eval_core(dataset, mode, scope)
%LOOPER_EVAL_CORE Shared evaluator for Atanas + Kato LOOPER runs.
%   dataset: "atanas" | "kato"
%   mode:    "stationarity" | "fidelity"
%   scope:   "single" | "all"

    if nargin < 3
        error('looper_eval_core requires dataset, mode, and scope.');
    end
    dataset = lower(string(dataset));
    mode = lower(string(mode));
    scope = lower(string(scope));
    if dataset ~= "atanas" && dataset ~= "kato"
        error('dataset must be "atanas" or "kato".');
    end
    if mode ~= "stationarity" && mode ~= "fidelity"
        error('mode must be "stationarity" or "fidelity".');
    end
    if scope ~= "single" && scope ~= "all"
        error('scope must be "single" or "all".');
    end

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir); % run_looper_diagnostics + helpers
    if dataset == "atanas"
        addpath(fullfile(rootDir, 'atanas-data'));
    else
        addpath(fullfile(rootDir, 'kato_2015'));
        addpath(fullfile(rootDir, 'kato_looper'));
    end

    % Ensure figures open as docked tabs.
    set(0, 'DefaultFigureWindowStyle', 'docked');
    set(0, 'DefaultFigureVisible', 'on');

    outDir = results_dir_(rootDir, dataset, scope, mode);
    if mode == "stationarity"
        plotDir = fullfile(outDir, 'plots');
        if ~exist(plotDir, 'dir')
            mkdir(plotDir);
        end
        projDir = fullfile(outDir, 'projections');
        if ~exist(projDir, 'dir')
            mkdir(projDir);
        end
    end

    if scope == "single"
        matPath = fullfile(outDir, sprintf('%s_%s_%s.mat', dataset, scope, mode));
        if ~exist(matPath, 'file')
            error('No %s %s results found. Run %s_%s_%s_run first.', ...
                dataset, mode, dataset, scope, mode);
        end

        S = load(matPath);
        if ~isfield(S, 'saveData')
            error('Expected saveData in %s', matPath);
        end
        if ~isfield(S, 'worm')
            error('Expected worm metadata in %s', matPath);
        end

        worm = S.worm;
        saveData = S.saveData;
        safeUid = regexprep(worm.uid, '[^A-Za-z0-9_-]', '_');
        tag = sprintf('worm_%s', safeUid);

        dt = worm.dt_sec;
        if isempty(dt) || ~isfinite(dt)
            dt = 1;
        end

        diagDir = fullfile(outDir, 'diagnostics');
        if ~exist(diagDir, 'dir')
            mkdir(diagDir);
        end

        if ~has_required_fields_(saveData)
            fprintf('Skipping %s (missing LOOPER fields).\n', tag);
            meta = meta_from_worm_(worm, dt);
            meta.status = 'failed_missing_fields';
            summary = looper_helpers.base_summary(struct(), struct(), meta);
        elseif mode == "stationarity"
            [preIdx, postIdx, t_on] = stationarity_split_(dataset, worm, saveData);
            diag = run_looper_diagnostics(saveData, struct('tag', tag, 'outDir', diagDir, ...
                'dt_sec', dt, 'preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));

            meta = meta_from_worm_(worm, dt);
            evalOut = looper_helpers.eval_stationarity(worm.RawData, saveData, ...
                preIdx, postIdx, dt, meta, diag);
            looper_helpers.save_stationarity_artifacts(evalOut, worm, plotDir, projDir, ...
                struct('tag', tag, 'makeLoopPlot', false));

            summary = evalOut.summary;
        else
            diag = run_looper_diagnostics(saveData, struct('tag', tag, 'outDir', diagDir, 'dt_sec', dt));
            meta = meta_from_worm_(worm, dt);
            summary = looper_helpers.base_summary(saveData, diag, meta);
        end

        writetable(struct2table(summary), fullfile(outDir, 'summary.csv'));
        fprintf('Diagnostics saved in: %s\n', diagDir);
        return;
    end

    if ~exist(outDir, 'dir')
        error('No results folder found: %s. Run %s_%s_%s_run first.', outDir, dataset, scope, mode);
    end
    files = dir(fullfile(outDir, 'worm_*.mat'));
    if isempty(files)
        error('No %s %s results found in %s. Run %s_%s_%s_run first.', ...
            dataset, mode, outDir, dataset, scope, mode);
    end

    diagRoot = fullfile(outDir, 'diagnostics');
    if ~exist(diagRoot, 'dir')
        mkdir(diagRoot);
    end

    summaries = struct([]);
    for i = 1:numel(files)
        matPath = fullfile(files(i).folder, files(i).name);
        S = load(matPath);
        if ~isfield(S, 'saveData') || ~isfield(S, 'worm')
            fprintf('Skipping %s (missing worm/saveData).\n', files(i).name);
            continue;
        end

        worm = S.worm;
        saveData = S.saveData;
        safeUid = regexprep(worm.uid, '[^A-Za-z0-9_-]', '_');
        tag = sprintf('worm_%s', safeUid);

        dt = worm.dt_sec;
        if isempty(dt) || ~isfinite(dt)
            dt = 1;
        end

        diagDir = fullfile(diagRoot, tag);
        if ~has_required_fields_(saveData)
            fprintf('Skipping %s (missing LOOPER fields).\n', tag);
            meta = meta_from_worm_(worm, dt);
            meta.status = 'failed_missing_fields';
            summary = looper_helpers.base_summary(struct(), struct(), meta);
        elseif mode == "stationarity"
            [preIdx, postIdx, t_on] = stationarity_split_(dataset, worm, saveData);
            diag = run_looper_diagnostics(saveData, struct('tag', tag, 'outDir', diagDir, ...
                'dt_sec', dt, 'preIdxRaw', preIdx, 'postIdxRaw', postIdx, 'tOnRaw', t_on));

            meta = meta_from_worm_(worm, dt);
            evalOut = looper_helpers.eval_stationarity(worm.RawData, saveData, ...
                preIdx, postIdx, dt, meta, diag);
            looper_helpers.save_stationarity_artifacts(evalOut, worm, plotDir, projDir, ...
                struct('tag', tag, 'makeLoopPlot', true));

            summary = evalOut.summary;
        else
            diag = run_looper_diagnostics(saveData, struct('tag', tag, 'outDir', diagDir, 'dt_sec', dt));
            meta = meta_from_worm_(worm, dt);
            summary = looper_helpers.base_summary(saveData, diag, meta);
        end

        if isempty(summaries)
            summaries = summary;
        else
            summaries(end+1) = summary; %#ok<SAGROW>
        end
    end

    if ~isempty(summaries)
        writetable(struct2table(summaries), fullfile(outDir, 'summary.csv'));
    end
    fprintf('Diagnostics saved in: %s\n', diagRoot);
end

function outDir = results_dir_(rootDir, dataset, scope, mode)
    outDir = fullfile(rootDir, 'results', sprintf('%s_%s', dataset, scope), mode);
end

function meta = meta_from_worm_(worm, dt)
    meta = struct('uid', worm.uid, 'n_neurons', worm.num_neurons, 'T', worm.T, 'dt_sec', dt);
end

function [preIdx, postIdx, t_on] = stationarity_split_(dataset, worm, saveData)
    if ~isfield(saveData, 'TrainSplit') || ~isfield(saveData.TrainSplit, 'enabled') ...
            || ~saveData.TrainSplit.enabled
        error('Expected TrainSplit for stationarity evaluation.');
    end
    preIdx = saveData.TrainSplit.preIdx;
    postIdx = saveData.TrainSplit.postIdx;
    t_on = postIdx(1);
end

function ok = has_required_fields_(saveData)
    required = {'BestStateMap', 'BestLoopAssignments', 'FinalStream'};
    ok = isstruct(saveData);
    for i = 1:numel(required)
        if ~ok || ~isfield(saveData, required{i}) || isempty(saveData.(required{i}))
            ok = false;
            return;
        end
    end
end
